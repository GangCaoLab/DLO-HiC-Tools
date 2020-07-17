from os.path import join, split, splitext, dirname, basename
import os.path as osp
from copy import copy
import logging
import multiprocessing as mp
import subprocess
from collections import OrderedDict
import pysam
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import re
import os
import glob

import click
from pysam import Samfile

from dlo_hic.utils.wrap.bwa import BWA
from dlo_hic.utils.parse_text import Bedpe
from dlo_hic.utils.fastqio import Fastq, write_fastq
from dlo_hic.utils.stream import beds2bedpe, upper_triangle, write_to_file


log = logging.getLogger(__name__)


def process_sam(bwa, output, next_iter_fq=None, mapq_thresh=20, sam=None):
    """
    Process alignment results.

    Count alignment items:
        1. unique mapped (mapq >= mapq_thresh)
        2. unmapped (flag & 4 != 0)
        2. other (multiple mapped)

    Parameters
    ----------
    bwa : BWA
        bwa command wraper
    output : str
        Path to output BED6 file.
    next_iter_fq : str, optional
        Path to first round iteration fastq.
    mapq_thresh : int
        mapq threshold.
    sam : str, optional
        If keep sam, specify path to BAM file.

    Return
    ------
    counts : dict
        Number of three kinds of alignment result of input pair.
    """
    counts = {
        "unique": 0,
        "unmapped": 0,
        "other": 0,
    }
    sam_iter = pysam.AlignmentFile(bwa.samse(stdout=True), mode='rb')

    if next_iter_fq:
        # create iter records and fastq file for iterative mapping
        fq_fh = open(next_iter_fq, 'w')
    if sam:
        sam_fh_u = pysam.AlignmentFile(re.sub(".sam$", ".uniq.sam", sam), 'w', template=sam_iter)
        sam_fh_n = pysam.AlignmentFile(re.sub(".sam$", ".un.sam"  , sam), 'w', template=sam_iter)
        sam_fh_m = pysam.AlignmentFile(re.sub(".sam$", ".mul.sam" , sam), 'w', template=sam_iter)

    with open(output, 'w') as out:
        for read in sam_iter:

            s_type = sam_read_type(read, mapq_thresh)
            counts[s_type] += 1

            if sam:
                if s_type == 'unique':
                    sam_fh_u.write(read)
                elif s_type == 'unmapped':
                    sam_fh_n.write(read)
                else:
                    sam_fh_m.write(read)

            if s_type == 'unique':  # unique alignment
                bed = sam_to_bed(read)
                out.write(bed + '\n')
            elif next_iter_fq:  # write fastq files for iterative mapping
                if s_type == 'unmapped':
                    fq = sam_to_fq(read)
                    write_fastq(fq, fq_fh)

    if next_iter_fq:
        fq_fh.close()
    if sam:
        for fh in sam_fh_u, sam_fh_n, sam_fh_m:
            fh.close()

    return counts


def sam_read_type(read, mapq_thresh):
    """ Determine sam read type """
    if (read.mapq >= mapq_thresh) & (read.flag & 4 == 0):  # unique map
        return 'unique'
    elif read.flag & 4 != 0:
        return 'unmapped'
    else:
        return 'other'


def sam_to_fq(read):
    fq = Fastq(read.qname, read.seq, read.qual)
    return fq


def sam_to_bed(read):
    name = read.qname
    chr_, start, end = read.reference_name, read.reference_start, read.reference_end
    strand = '-' if read.is_reverse else '+'
    score = read.mapq
    bed_line = "\t".join([chr_, str(start), str(end), name, str(score), strand])
    return bed_line


def log_alignment_info(pet1_count, pet2_count, log_file=None):
    def count_msg(counts):
        total = sum(counts.values())
        u = counts['unique']
        m = counts['unmapped']
        o = counts['other']
        u_r = 0 if total == 0 else u/total
        m_r = 0 if total == 0 else m/total
        o_r = 0 if total == 0 else o/total
        msg = "unique-mapped\t{}\tpercent\t{:.2%}".format(u, u_r) + "\t"
        msg += "unmapped-mapped\t{}\tpercent\t{:.2%}".format(m, m_r) + "\t"
        msg += "multiple-mapped\t{}\tpercent\t{:.2%}".format(o, o_r)
        return (u, m, o), msg
    c1, msg1 = count_msg(pet1_count)
    c2, msg2 = count_msg(pet2_count)
    log.info("PET1 alignment results:\t" + msg1)
    log.info("PET2 alignment results:\t" + msg2)
    if log_file:
        with open(log_file, 'w') as fo:
            fo.write("# PET1\n")
            fo.write("unique-mapped\t{}\n".format(c1[0]))
            fo.write("unmapped-mapped\t{}\n".format(c1[1]))
            fo.write("multiple-mapped\t{}\n".format(c1[2]))
            fo.write("\n")
            fo.write("# PET2\n")
            fo.write("unique-mapped\t{}\n".format(c2[0]))
            fo.write("unmapped-mapped\t{}\n".format(c2[1]))
            fo.write("multiple-mapped\t{}\n".format(c2[2]))


def log_paired(counts):
    total = counts['total']
    ratio = 0 if total == 0 else counts['paired'] / total
    log.info("Unique paired reads: {}\tpercent: {:.2%}".format(counts['paired'], ratio))


@click.command(name="build_bedpe")
@click.argument("input1", nargs=1)
@click.argument("input2", nargs=1)
@click.argument("output", nargs=1)
@click.option("--ncpu", "-p",
    type=int,
    default=mp.cpu_count(),
    help="how many threads used to run bwa")
@click.option("--bwa-index",
    required=True,
    help="The bwa index prefix.")
@click.option("--mapq",
    default=1,
    help="the mapq threshold used to filter mapped records. default 1(for fetch unique mapping)")
@click.option("--iteration",
    default=1,
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
@click.option("--keep-sam", 
    is_flag=True,
    help="Keep SAM file or not.")
def _main(input1, input2, output, ncpu, bwa_index, mapq, iteration, log_file, bwa_log_file, keep_iter_files, keep_sam):
    """
    Build bedpe file from fastq files.
    
    \b
    Arguments
    ---------
    input1 : str
        Path to PET1 fastq file. 
    input2 : str
        Path to PET2 fastq file.
    output : str
        Output prefix.
    """
    log.info("Build BEDPE from %s %s"%(input1, input2))

    out_prefix = output

    iter_id = 0
    bwa1 = BWA(bwa_index, algorithm='aln', log_file=bwa_log_file+'.pet1') 
    bwa2 = BWA(bwa_index, algorithm='aln', log_file=bwa_log_file+'.pet2') 
    counts = {}
    pool = ProcessPoolExecutor(max_workers=2)

    log.info("Begin iterative mapping.")
    assert iteration >= 1, "iteration at least 1, {} got.".format(iteration)
    while iter_id < iteration:

        if iter_id == 0:
            fq1 = input1
            fq2 = input2
        else:
            fq1 = next_fq1
            fq2 = next_fq2

        log.info("-"*30)
        log.info("Iteration [{}]".format(iter_id))
        log.info("Do alignment using 'bwa aln', with {} threads".format(ncpu))

        # modify max diff parameter in each iteration
        bwa1.run(fq1, out_prefix+".pet1", thread=ncpu, max_diff=iter_id)
        bwa2.run(fq2, out_prefix+".pet2", thread=ncpu, max_diff=iter_id)

        log.info("'bwa aln' done, process the output")

        if iter_id >= iteration - 1:
            next_fq1 = next_fq2 = None
        else:
            next_fq1 = out_prefix + '.iter.' + str(iter_id+1) + '.pet1.fq'
            next_fq2 = out_prefix + '.iter.' + str(iter_id+1) + '.pet2.fq'

        out_bed1 = out_prefix + '.iter.{}.pet1.bed'.format(iter_id)
        out_bed2 = out_prefix + '.iter.{}.pet2.bed'.format(iter_id)

        if keep_sam:
            sam_1 = out_prefix + '.iter.{}.pet1.sam'.format(iter_id)
            sam_2 = out_prefix + '.iter.{}.pet2.sam'.format(iter_id)
        else:
            sam_1 = sam_2 = None

        # process pet1 and 2 in parallel
        map_ = pool.map 
#        map_ = map
        args = [bwa1, bwa2], [out_bed1, out_bed2], [next_fq1, next_fq2], repeat(mapq), [sam_1, sam_2]
        counts1, counts2 = map_(process_sam, *args)

        counts_iter = {
            'pet1': counts1,
            'pet2': counts2,
        }
        counts['iter{}'.format(iter_id)] = counts_iter
        log_alignment_info(counts_iter['pet1'], counts_iter['pet2'])

        iter_id += 1

        if ((counts_iter['pet1']['unique']    == 0) and (counts_iter['pet2']['unique']   == 0)) or \
           ((counts_iter['pet1']['unmapped']  == 0) and (counts_iter['pet2']['unmapped'] == 0)):
            # break loop when no new unique or no unmatched found in this iteration
            break

    log.info("-"*30)
    log.info("Iterative mapping stopped.")

    # merge files
    outdir = dirname(osp.abspath(out_prefix))
    files_in_outdir = os.listdir(outdir)

    file_collections = {
        (tp, pet): glob.glob(join(out_prefix+".iter.*."+pet+"."+tp))
        for tp in ['bed'] + (['uniq.sam', 'un.sam', 'mul.sam'] if keep_sam else [])
        for pet in ['pet1', 'pet2']
    }
    for (tp, pet), files in file_collections.items():
        res = out_prefix + '.' + pet + '.' + tp
        if len(files) > 1:
            if 'sam' in tp:
                cmd = "cat " + " ".join(files) + " | grep -v @  > " + res
            else:
                cmd = "cat " + " ".join(files) + " > " + res
            subprocess.check_call(cmd, shell=True)
        else:
            subprocess.check_call(['mv', files[0], res])

    bed1 = out_prefix + '.pet1.bed'
    bed2 = out_prefix + '.pet2.bed'
    bed1_sorted = re.sub(".bed$", ".sorted.bed", bed1)
    bed2_sorted = re.sub(".bed$", ".sorted.bed", bed2)
    # sort BED files by seq id
    subprocess.check_call("sort --parallel={} -k4,4 {} > {}".format(ncpu, bed1, bed1_sorted), shell=True)
    subprocess.check_call("sort --parallel={} -k4,4 {} > {}".format(ncpu, bed2, bed2_sorted), shell=True)
    # clean un sorted bed
    subprocess.check_call(['rm', bed1])
    subprocess.check_call(['rm', bed2])
    # join pet1 and 2
    bedpe = out_prefix + '.uniq.bedpe'
    msg = "join {} {} to bedpe: {}".format(bed1_sorted, bed2_sorted, bedpe)
    log.info(msg)
    line_itr = beds2bedpe(bed1_sorted, bed2_sorted)
    line_itr = upper_triangle(line_itr, 'bedpe')
    line_cnt = write_to_file(line_itr, bedpe)

    counts['final'] = {
        'paired': line_cnt,
        'pet1': {},
        'pet2': {},
    }

    # aggregate the counts of each iteration
    counts_iters = {
        int(k.replace('iter','')):d
            for k,d in counts.items() if k.startswith('iter')}
    for p in 'pet1', 'pet2':
        counts['final'][p]['unique']   = sum([d[p]['unique']   for d in counts_iters.values()])
        counts['final'][p]['other']    = sum([d[p]['other']    for d in counts_iters.values()])
        counts['final'][p]['unmapped'] = counts_iters[sorted(counts_iters.keys())[-1]][p]['unmapped']

    log_alignment_info(counts['final']['pet1'], counts['final']['pet2'], log_file)

    total = sum(counts['final']['pet1'].values())
    ratio = 0 if total == 0 else line_cnt/total
    log.info("Unique paired reads:\t{}\tpercent\t{:.2%}".format(line_cnt, ratio))
    with open(log_file, "a") as f:
        f.write("\n")
        f.write("# Unique paired\n")
        f.write("total\t{}\n".format(line_cnt))

    # clean intermediate files
    subprocess.check_call(["rm", bed1_sorted, bed2_sorted])
    subprocess.check_call("rm " + out_prefix + ".*.sai", shell=True)
    if not keep_iter_files:
        if iteration > 1:
            log.info("Clean all iteration files.")
            subprocess.check_call("rm " + out_prefix + "*.iter.*", shell=True)

    log.info("Build BEDPE done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")