from os.path import join, split, splitext, dirname
from copy import copy
import logging
import multiprocessing as mp
import subprocess
from collections import OrderedDict

import click
from pysam import Samfile

from dlo_hic.utils.wrap.bwa import BWA
from dlo_hic.utils.parse_text import Bedpe
from dlo_hic.utils.fastqio import Fastq, write_fastq


log = logging.getLogger(__name__)


class IterRec(object):
    """ Abstraction of an iteration intermediate record. """
    def __init__(self, iter_id,
                 seqid, seq1, seq2, align_tp1, align_tp2,
                 chr1, start1, end1, chr2, start2, end2, strand1, strand2):
        self.iter_id = iter_id
        self.seqid = seqid
        self.seq1, self.seq2 = seq1, seq2
        self.align_tp1, self.align_tp2 = align_tp1, align_tp2
        self.chr1, self.start1, self.end1 = chr1, start1, end1
        self.chr2, self.start2, self.end2 = chr2, start2, end2
        if chr1 and start1 and end1:
            self.start1 = int(self.start1)
            self.end1 = int(self.end1)
        if chr2 and start2 and end2:
            self.start2 = int(self.start2)
            self.end2 = int(self.end2)
        self.strand1, self.strand2 = strand1, strand2

    def __str__(self):
        s1 = str(self.start1) if self.start1 else ""
        e1 = str(self.end1) if self.end1 else ""
        s2 = str(self.start2) if self.start2 else ""
        e2 = str(self.end2) if self.end2 else ""
        fields = [self.seqid, self.seq1, self.seq2, self.align_tp1, self.align_tp2,
                  self.chr1 or "", s1, e1, self.chr2 or "", s2, e2,
                  self.strand1 or "", self.strand2 or ""]
        return "\t".join(fields)

    @classmethod
    def from_line(cls, line, iter_id):
        """ Construct from a string. """
        items = line.split("\t")
        rec = cls(iter_id, *items)
        return rec

    @staticmethod
    def strand_of_sam_rec(rec):
        if rec.reference_name:
            strand = '-' if rec.is_reverse else '+'
        else:
            strand = None
        return strand

    @classmethod
    def from_sam_pair(cls, rec1, rec2, iter_id, mapq_thresh):
        """ Construct from sam record pair. """
        seqid = rec1.query_name
        seq1, seq2 = rec1.seq, rec2.seq
        align_tp1 = sam_read_type(rec1, mapq_thresh)
        align_tp2 = sam_read_type(rec2, mapq_thresh)
        chr1 = rec1.reference_name
        start1, end1 = rec1.reference_start, rec1.reference_end
        chr2 = rec2.reference_name
        start2, end2 = rec2.reference_start, rec2.reference_end
        strand1 = IterRec.strand_of_sam_rec(rec1)
        strand2 = IterRec.strand_of_sam_rec(rec2)
        rec = cls(iter_id, seqid, seq1, seq2, align_tp1, align_tp2,
                  chr1, start1, end1, chr2, start2, end2, strand1, strand2)
        return rec

    @staticmethod
    def get_subseqs(seq, iterid):
        """ Get all sub-sequences """
        sub_seqs = []
        sub_len = len(seq) - iterid
        for start in range(0, iterid+1):
            end = start + sub_len
            sub = seq[start:end]
            sub_seqs.append(sub)
        return sub_seqs

    def to_fq_recs(self, quality_chr='A'):
        """ Get all sub-sequences's fastq records. """
        fq_recs = []
        for seq, aln_tp in (self.seq1, self.align_tp1), (self.seq2, self.align_tp2):
            pet_label = 'pet1' if seq is self.seq1 else 'pet2'
            if aln_tp != 'unique':  # only convert non-unique mapped read
                for idx, subseq in enumerate(IterRec.get_subseqs(seq, self.iter_id)):
                    seq_id = self.seqid + '_{}_'.format(pet_label) + str(idx)
                    qua = quality_chr * len(subseq)
                    fq = Fastq(seq_id, subseq, qua)
                    fq_recs.append(fq)
        return fq_recs

    @property
    def is_unique_paired(self):
        return self.align_tp2 == self.align_tp1 == "unique"

    def to_bedpe(self):
        """ Convert to BEDPE record. """
        bpe = Bedpe(self.chr1, self.start1, self.end1, self.chr2, self.start2, self.end2,
                    self.seqid, 0, self.strand1, self.strand2)
        return bpe


def process_sam_pair(input1, input2, output, iter1_fq=None, iter_db=None, mapq_thresh=20):
    """
    Process alignment result in first iteration.
    Merge the unique pair to bedpe file, and store unpaired data to an intermediate file.

    Count alignment items:
        1. unique mapped (mapq >= mapq_thresh & not have 'XA' tag)
        2. multiple mapped (have 'XA' tag)
        3. other (unmatched or too many multiple match)

    Parameters
    ----------
    input1 : str
        Path to input1 SAM/BAM file.
    input2 : str
        Path to input2 SAM/BAM file.
    output : str
        Path to output BEDPE file.
    iter1_fq : str, optional
        Path to first round iteration fastq.
    iter_db : dict, optional
        iter dict
    mapq_thresh : int
        mapq threshold.

    Return
    ------
    counts : dict
        Number of three kinds of alignment result of input pair.
    """
    prefix1, fmt1 = splitext(input1)
    mode_r1 = "rb" if fmt1 == ".bam" else "r"
    prefix2, fmt2 = splitext(input2)
    mode_r2 = "rb" if fmt2 == ".bam" else "r"

    count_items = {
        "unique": 0,
        "multiple": 0,
        "other": 0,
    }
    counts = {
        'input1': copy(count_items),
        'input2': copy(count_items),
        'paired': 0,
        'total': 0,
    }

    if iter1_fq:
        # create iter records and fastq file for iterative mapping
        fq_fh = open(iter1_fq, 'w')

    with Samfile(input1, mode_r1) as sam1, Samfile(input2, mode_r2) as sam2,\
            open(output, 'w') as out:
        iter1 = sam1.fetch(until_eof=True)
        iter2 = sam2.fetch(until_eof=True)
        while True:
            try:
                read1 = next(iter1)
                read2 = next(iter2)
            except StopIteration:
                break

            r1_type = sam_read_type(read1, mapq_thresh)
            r2_type = sam_read_type(read2, mapq_thresh)

            counts['input1'][r1_type] += 1
            counts['input2'][r2_type] += 1

            is_paired = (r1_type == r2_type == 'unique')

            if is_paired:  # unique pair, write to BEDPE
                counts['paired'] += 1
                bpe = sam_pair_to_bedpe(read1, read2)
                out.write(str(bpe)+'\n')
            else:
                if iter1_fq:  # write files for iterative mapping
                    iter_rec = IterRec.from_sam_pair(read1, read2, 1, mapq_thresh)
                    iter_db[iter_rec.seqid] = iter_rec
                    for fq in iter_rec.to_fq_recs():
                        write_fastq(fq, fq_fh)

            counts['total'] += 1

    if iter1_fq:
        fq_fh.close()

    return counts


def sam_read_type(read, mapq_thresh):
    """ Determine sam read type """
    def has_xa(read):
        return True if len(read.tags) > 0 and read.tags[-1][0] == 'XA' else False

    xa = has_xa(read)
    if (read.mapq >= mapq_thresh) and (not xa):  # unique map
        return 'unique'
    elif xa:
        return 'multiple'
    else:
        return 'other'


def sam_pair_to_bedpe(read1, read2):
    """ merge samme reads pair to BEDPE record. """
    chr1, start1, end1 = read1.reference_name, read1.reference_start, read1.reference_end
    chr2, start2, end2 = read2.reference_name, read2.reference_start, read2.reference_end
    strand1 = '-' if read1.is_reverse else '+'
    strand2 = '-' if read2.is_reverse else '+'
    bpe = Bedpe(chr1, start1, end1, chr2, start2, end2, read1.query_name, 0, strand1, strand2)
    bpe.to_upper_trangle()
    return bpe


def iterative_mapping(bwa_index, iter_db, iter1_fq, bedpe, n_iters, counts, threads, mapq_thresh):
    """
    Take the sub-sequences of non-unique mapped reads,
    alignment again, increase the reads use ratio.

    Parameters
    ----------
    bwa_index : str
        Prefix of BWA index files.
    iter_db : dict
        dict for store unpaired records.
    iter1_fq : str
        Path to Fastq file in first iteration.
    bedpe : str
        Path to output BEDPE file.
    n_iters : int
        Number of iteration.
    counts : dict
        Counter.
    threads : int
        Number of threads for run BWA.
    mapq_thresh : int
        mapq threshold.
    """
    bpe_fh = open(bedpe, 'a')
    counts['iters'] = {}

    def process_sam_records_belong_to_same_read(iter_id, records, seq_id, pet_label):
        unique_recs = []
        for rec in records:
            if sam_read_type(rec, mapq_thresh) == 'unique':
                unique_recs.append(rec)
        if len(unique_recs) == 0:
            return
        elif len(unique_recs) > 1:
            if seq_id in iter_db:
                del iter_db[seq_id]
        else:  # only process reads, which only one sub-seq unique mapped
            u_rec = unique_recs[0]
            try:
                iter_rec = iter_db[seq_id]
            except KeyError:  # seqid already clean out
                return
            if pet_label == 'pet1':
                iter_rec.align_tp1 = "unique"
                iter_rec.chr1 = u_rec.reference_name
                iter_rec.start1 = u_rec.reference_start
                iter_rec.end1 = u_rec.reference_end
                iter_rec.strand1 = '-' if u_rec.is_reverse else '+'
            else:
                iter_rec.align_tp2 = "unique"
                iter_rec.chr2 = u_rec.reference_name
                iter_rec.start2 = u_rec.reference_start
                iter_rec.end2 = u_rec.reference_end
                iter_rec.strand2 = '-' if u_rec.is_reverse else '+'
            iter_db[seq_id] = iter_rec  # update db
            if iter_rec.is_unique_paired:  # is paired, delete record, write pair to bedpe
                counts['iters'][iter_id]['paired'] += 1
                counts['paired'] += 1
                bpe_fh.write(str(iter_rec.to_bedpe()) + "\n")
                del iter_db[seq_id]

    def reproduce_fastq(fq_fh, iter_id):
        for seq_id, rec in iter_db.items():
            for fq in rec.to_fq_recs():
                write_fastq(fq, fq_fh)

    for iter_id in range(1, n_iters+1):
        counts['iters'][iter_id] = {'paired': 0}
        log.info("-"*30)
        log.info("Iterative mapping, <round {}>".format(iter_id))
        iter_fq = bedpe + '.iter.{}.fq'.format(iter_id)

        log.info("Do alignment using 'bwa aln', with {} threads".format(threads))
        bwa_log_file = bedpe + '.iter.' + str(iter_id) + '.bwa.log'
        outdir = dirname(bedpe)
        pre = join(outdir, splitext(split(iter_fq)[1])[0])  # prefix of input fastq file
        bwa = BWA(bwa_index, algorithm='aln', log_file=bwa_log_file)
        bwa.run(iter_fq, pre, thread=threads, max_diff=0, bam=True)
        bam = pre + '.bam'

        # process bam
        log.info("Process BAM file")
        with Samfile(bam, 'rb') as sam:
            sam_iter = sam.fetch(until_eof=True)
            old = None
            same_read_recs = []
            for rec in sam_iter:
                rec_id = rec.query_name
                seq_id, pet_label, r_idx = rec_id.split("_")
                if (seq_id, pet_label) == old:  # read id same to previous one
                    same_read_recs.append(rec)
                elif old is not None:  # new read id, process previous records
                    process_sam_records_belong_to_same_read(iter_id, same_read_recs, old[0], old[1])
                    same_read_recs = [rec]
                old = (seq_id, pet_label)
            process_sam_records_belong_to_same_read(iter_id, same_read_recs, old[0], old[1])
        
        new_uniq_pairs = counts['iters'][iter_id]['paired']
        log.info("{} new unique pairs found in <round {}>".format(new_uniq_pairs, iter_id))
        log_paired(counts)

        if new_uniq_pairs == 0:
            log.info("No new unique pairs found, stop iteration.")
            break

        if iter_id+1 <= n_iters:
            new_fq = bedpe + '.iter.{}.fq'.format(iter_id+1)
            with open(new_fq, 'w') as fq_fh:
                log.info("Reproduce FASTQ: {}".format(new_fq))
                reproduce_fastq(fq_fh, iter_id+1)

    bpe_fh.close()


def log_alignment_info(pet1_count, pet2_count, log_file=None):
    def count_msg(counts):
        total = sum(counts.values())
        u = counts['unique']
        m = counts['multiple']
        o = counts['other']
        u_r = 0 if total == 0 else u/total
        m_r = 0 if total == 0 else m/total
        o_r = 0 if total == 0 else o/total
        msg = "unique-mapped\t{}\tpercent\t{:.2%}".format(u, u_r) + "\t"
        msg += "multiple-mapped\t{}\tpercent\t{:.2%}".format(m, m_r) + "\t"
        msg += "other\t{}\tpercent\t{:.2%}".format(o, o_r)
        return (u, m, o), msg
    c1, msg1 = count_msg(pet1_count)
    c2, msg2 = count_msg(pet2_count)
    log.info("PET1 alignment results:\t" + msg1)
    log.info("PET2 alignment results:\t" + msg2)
    if log_file:
        with open(log_file, 'w') as fo:
            fo.write("# PET1\n")
            fo.write("unique-mapped\t{}\n".format(c1[0]))
            fo.write("multiple-mapped\t{}\n".format(c1[1]))
            fo.write("other\t{}\n".format(c1[2]))
            fo.write("\n")
            fo.write("# PET2\n")
            fo.write("unique-mapped\t{}\n".format(c2[0]))
            fo.write("multiple-mapped\t{}\n".format(c2[1]))
            fo.write("other\t{}\n".format(c2[2]))


def log_paired(counts):
    total = counts['total']
    ratio = 0 if total == 0 else counts['paired'] / total
    log.info("Unique paired reads: {}\tpercent: {:.2%}".format(counts['paired'], ratio))


@click.command(name="build_bedpe")
@click.option("--file-format", "-f",
    default="fastq",
    type=click.Choice(['fastq', 'bam', 'sam']))
@click.argument("input1", nargs=1)
@click.argument("input2", nargs=1)
@click.argument("bedpe", nargs=1)
@click.option("--ncpu", "-p",
    type=int,
    default=mp.cpu_count(),
    help="how many threads used to run bwa")
@click.option("--bwa-index", help="The bwa index prefix. If the input format is fastq, it's required.")
@click.option("--mapq",
    default=1,
    help="the mapq threshold used to filter mapped records. default 1(for fetch unique mapping)")
@click.option("--iterative",
    default=0,
    help="Iterative alignment, specify the number of iterative loop.")
@click.option("--log-file",
    default="build_bedpe.log",
    help="Sperate log file for storage some count information.")
@click.option("--bwa-log-file",
    default="bwa.log",
    help="Separate log file for storage bwa output")
@click.option("--keep-iter-files",
    is_flag=True,
    help="Don't delete intermediate iterations files")
def _main(file_format, input1, input2, bedpe, ncpu, bwa_index, mapq, iterative, log_file, bwa_log_file, keep_iter_files):
    """ Build bedpe file from fastq or sam/bam file. """
    log.info("Build BEDPE from %s %s"%(input1, input2))

    outdir = dirname(bedpe)

    pre_1 = join(outdir, splitext(split(input1)[1])[0]) # prefix of input fastq file
    pre_2 = join(outdir, splitext(split(input2)[1])[0])
    if file_format == 'fastq':
        log.info("Do alignment using 'bwa aln', with {} threads".format(ncpu))
        assert bwa_index is not None, 'bwa_index is required when input fastq file.'
        bwa = BWA(bwa_index, algorithm='aln', log_file=bwa_log_file) # write bwa log to separate file
        bwa.run(input1, pre_1, thread=ncpu, max_diff=0, bam=True)
        bwa.run(input2, pre_2, thread=ncpu, max_diff=0, bam=True)
        bam1 = pre_1 + '.bam'
        bam2 = pre_2 + '.bam'
    else:
        # build from bam file directly
        bam1 = input1
        bam2 = input2

    log.info("Process BAM files.")
    if iterative > 0:
        log.info("Initialize for iterative mapping.")
        assert bwa_index is not None, 'bwa_index is required, when do iterative mapping.'
        iter1_fq = bedpe + '.iter.1.fq'
        iter_db = OrderedDict()
    else:
        iter1_fq = iter_db = None

    counts = process_sam_pair(bam1, bam2, bedpe, iter1_fq, iter_db, mapq_thresh=mapq)
    log_alignment_info(counts['input1'], counts['input2'], log_file)
    log_paired(counts)

    if iterative > 0:
        log.info("Begin iterative mapping")
        iterative_mapping(bwa_index, iter_db, iter1_fq, bedpe, iterative, counts, ncpu, mapq)

        # clear intermedia files
        if not keep_iter_files:
            subprocess.check_call("rm {}.iter.*".format(bedpe), shell=True)

    with open(log_file, "a") as f:
        f.write("\n")
        f.write("# Unique paired\n")
        f.write("total\t{}\n".format(counts['paired']))

    log.info("Build BEDPE done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")