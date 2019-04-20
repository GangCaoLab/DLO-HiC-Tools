from os.path import join, split, splitext, dirname
from copy import copy
import logging
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import click
from pysam import Samfile

from dlo_hic.utils.wrap.bwa import BWA
from dlo_hic.utils.parse_text import Bedpe


log = logging.getLogger(__name__)


def process_sam_pair(input1, input2, output, intermed, mapq_thresh=20):
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
    intermed : str
        Path to intermediate file.
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

    def has_xa(read):
        return True if len(read.tags) > 0 and read.tags[-1][0] == 'XA' else False

    def read_type(read):
        xa = has_xa(read)
        if (read.mapq >= mapq_thresh) and (not xa):  # unique map
            return 'unique'
        elif xa:
            return 'multiple'
        else:
            return 'other'

    with Samfile(input1, mode_r1) as sam1, Samfile(input2, mode_r2) as sam2,\
            open(output, 'w') as out, open(intermed, 'w') as interm:
        iter1 = sam1.fetch(until_eof=True)
        iter2 = sam2.fetch(until_eof=True)
        while True:
            try:
                read1 = next(iter1)
                read2 = next(iter2)
            except StopIteration:
                break

            r1_type = read_type(read1)
            r2_type = read_type(read2)

            counts['input1'][r1_type] += 1
            counts['input2'][r2_type] += 1

            is_paired = (r1_type == r2_type == 'unique')

            if is_paired:  # unique pair, write to BEDPE
                counts['paired'] += 1
                bpe = sam_pair_to_bedpe(read1, read2)
                out.write(str(bpe)+'\n')
            else:
                #interm.write()
                pass

            counts['total'] += 1

    return counts


def sam_pair_to_bedpe(read1, read2):
    """ merge samme reads pair to BEDPE record. """
    chr1, start1, end1 = read1.reference_name, read1.reference_start, read1.reference_end
    chr2, start2, end2 = read2.reference_name, read2.reference_start, read2.reference_end
    strand1 = '-' if read1.is_reverse else '+'
    strand2 = '-' if read2.is_reverse else '+'
    bpe = Bedpe(chr1, start1, end1, chr2, start2, end2, read1.query_name, 0, strand1, strand2)
    bpe.to_upper_trangle()
    return bpe


def log_split_count(pet1_count, pet2_count, log_file):
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
    default=1,
    help="Iterative alignment, specify the number of iterative loop.")
@click.option("--log-file",
    default="build_bedpe.log",
    help="Sperate log file for storage some count information.")
@click.option("--bwa-log-file",
    default="bwa.log",
    help="Separate log file for storage bwa output")
def _main(file_format, input1, input2, bedpe, ncpu, bwa_index, mapq, iterative, log_file, bwa_log_file):
    """ Build bedpe file from fastq or sam/bam file. """
    log.info("Build bedpe from %s %s"%(input1, input2))

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

    log.info("merge beds to bedpe.")
    intermed_file = pre_1 + '.med'
    counts = process_sam_pair(bam1, bam2, bedpe, intermed_file, mapq_thresh=mapq)

    log_split_count(counts['input1'], counts['input2'], log_file)

    with open(log_file, "a") as f:
        f.write("\n")
        f.write("# Unique paired\n")
        f.write("total\t{}\n".format(counts['paired']))

    log.info("Build bedpe done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")