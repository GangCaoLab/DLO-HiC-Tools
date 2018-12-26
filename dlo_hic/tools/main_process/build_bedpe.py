from os.path import join, split, splitext, dirname
import logging
import subprocess
import multiprocessing as mp

import click

from dlo_hic.utils.wrap.bwa import BWA
from dlo_hic.utils.stream import beds2bedpe, upper_triangle, write_to_file


log = logging.getLogger(__name__)


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
@click.option("--bwa-index", required=True, help="The bwa index prefix.")
@click.option("--mapq",
    default=1,
    help="the mapq threshold used to filter mapped records. default 1(for fetch unique mapping)")
@click.option("--upper-tri/non-upper-tri",
    default=True,
    help="Convert bedpe to upper-triangular format or not, default True")
@click.option("--bwa-log-file",
    default="bwa.log",
    help="separate log file for storage bwa output: default 'bwa.log'")
def _main(file_format, input1, input2, bedpe, ncpu, bwa_index, mapq, upper_tri, bwa_log_file):
    """ Build bedpe file from fastq or sam/bam file. """
    log.info("Build bedpe from %s %s"%(input1, input2))

    outdir = dirname(bedpe)

    pre_1 = join(outdir, splitext(split(input1)[1])[0]) # prefix
    pre_2 = join(outdir, splitext(split(input2)[1])[0])
    if file_format == 'fastq':
        log.info("Do alignment using 'bwa aln', with {} threads".format(ncpu))
        assert bwa_index is not None
        bwa = BWA(bwa_index, log_file=bwa_log_file) # write bwa log to separate file
        bwa.run(input1, pre_1, thread=ncpu, mem=False)
        bwa.run(input2, pre_2, thread=ncpu, mem=False)
        bam1 = pre_1 + '.bam'
        bam2 = pre_2 + '.bam'
    else:
        # build from bam file directly
        bam1 = input1
        bam2 = input2

    log.info("filter bam file to fetch unique mapping reads.")
    p1 = subprocess.Popen("samtools view {} -b -q {} > {}.filtered.bam".format(bam1, mapq, pre_1), shell=True)
    p2 = subprocess.Popen("samtools view {} -b -q {} > {}.filtered.bam".format(bam2, mapq, pre_2), shell=True)
    p1.wait()
    p2.wait()

    log.info("convert bam to bed.")
    bed1_ = pre_1 + ".bed.tmp"
    bed2_ = pre_2 + ".bed.tmp"
    p1 = subprocess.Popen("bedtools bamtobed -i {} > {}".format(pre_1+'.filtered.bam', bed1_), shell=True)
    p2 = subprocess.Popen("bedtools bamtobed -i {} > {}".format(pre_2+'.filtered.bam', bed2_), shell=True)
    p1.wait()
    p2.wait()

    log.info("sort bed file")
    bed1 = pre_1 + ".bed"
    bed2 = pre_2 + ".bed"
    subprocess.check_call("sort --parallel={} -k4,4 {} > {}".format(ncpu, bed1_, bed1), shell=True)
    subprocess.check_call("rm {}".format(bed1_), shell=True)
    subprocess.check_call("sort --parallel={} -k4,4 {} > {}".format(ncpu, bed2_, bed2), shell=True)
    subprocess.check_call("rm {}".format(bed2_), shell=True)

    log.info("merge beds to bedpe.")
    msg = "join {} {} to bedpe: {}".format(bed1, bed2, bedpe)
    log.info(msg)
    line_itr = beds2bedpe(bed1, bed2)
    if upper_tri:
        log.info("Convert BEDPE to upper triangular form.")
        line_itr = upper_triangle(line_itr, 'bedpe')
    write_to_file(line_itr, bedpe)

    log.info("Build bedpe done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")