import os
import logging
import subprocess
import multiprocessing as mp

import click

from dlo_hic.utils.wrap.bwa import BWA


log = logging.getLogger(__name__)


def beds2bedpe(bed1, bed2, bedpe_filename):
    msg = "join {} {} to bedpe: {}".format(bed1, bed2, bedpe_filename)
    log.info(msg)
    cmd = "join -j 4 {} {}".format(bed1, bed2)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    with open(bedpe_filename, 'w') as f:
        #
        # format joined bed bedpe
        for line in p.stdout:
            line = line.decode("utf-8")            
            items = line.strip().split()
            name, chr_a, s_a, e_a, score_a, strand_a,\
                  chr_b, s_b, e_b, score_b, strand_b = items
            outitems = [chr_a, s_a, e_a,
                        chr_b, s_b, e_b,
                        name, '0', strand_a, strand_b]
            outline = "\t".join(outitems) + "\n"
            f.write(outline)


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
@click.option("--bwa-log-file",
    default="bwa.log",
    help="separate log file for storage bwa output: default 'bwa.log'")
def _main(file_format, input1, input2, bedpe, ncpu, bwa_index, mapq, bwa_log_file):
    """ Build bedpe file from fastq or sam/bam file. """
    log.info("Build bedpe from %s %s"%(input1, input2))

    pre_1 = os.path.splitext(input1)[0] # prefix
    pre_2 = os.path.splitext(input2)[0]
    if file_format == 'fastq':
        # alignment firstly if the input format is fastq
        assert bwa_index is not None
        bwa = BWA(bwa_index, log_file=bwa_log_file) # write bwa log to separate file
        bwa.run(input1, pre_1, thread=ncpu, mem=False)
        bwa.run(input2, pre_2, thread=ncpu, mem=False)
        bam1 = pre_1 + '.bam'
        bam2 = pre_2 + '.bam'
    else:
        bam1 = input1
        bam2 = input2
    log.info("filter bam file to fetch unique mapping reads.")
    subprocess.check_call("samtools view {} -b -q {} > {}.filtered.bam".format(bam1, mapq, pre_1), shell=True)
    subprocess.check_call("samtools view {} -b -q {} > {}.filtered.bam".format(bam2, mapq, pre_2), shell=True)
    log.info("convert bam to bed.")
    bed1 = pre_1 + ".bed"
    bed2 = pre_2 + ".bed"
    subprocess.check_call("bedtools bamtobed -i {} | sort -k4,4 > {}".format(pre_1+'.filtered.bam', bed1), shell=True)
    subprocess.check_call("bedtools bamtobed -i {} | sort -k4,4 > {}".format(pre_2+'.filtered.bam', bed2), shell=True)
    log.info("merge beds to bedpe.")
    beds2bedpe(bed1, bed2, bedpe)


main = _main.callback


if __name__ == "__main__":
    _main()