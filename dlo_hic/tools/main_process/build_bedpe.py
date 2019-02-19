from os.path import join, split, splitext, dirname
import logging
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

import click
from pysam import Samfile

from dlo_hic.utils.wrap.bwa import BWA
from dlo_hic.utils.stream import beds2bedpe, upper_triangle, write_to_file


log = logging.getLogger(__name__)


def split_sam(input, mapq_thresh=20):
    """
    split alignment results(BAM/SAM) to:
        1. unique mapped (mapq >= mapq_thresh & not have 'XA' tag)
        2. multiple mapped (have 'XA' tag)
        3. other (unmatched or too many multiple match)

    Parameters
    ----------
    input : str
        Path to input SAM/BAM file.
    mapq_thresh : int
        mapq threshold.

    Return
    ------
    counts : dict
        Number of three kinds of alignment result.
    """
    prefix, fmt = splitext(input)
    uniq_file = prefix + ".uniq" + fmt
    mul_file = prefix + ".mul" + fmt
    unm_file = prefix + ".unm" + fmt
    mode_r = "rb" if fmt == ".bam" else "r"
    mode_w = "wb" if fmt == ".bam" else "w"
    counts = {
        "unique-mapped": 0,
        "multiple-mapped": 0,
        "other": 0,
    }
    with Samfile(input, mode_r) as sam:
        with Samfile(uniq_file, mode_w, template=sam) as uniq, \
             Samfile(mul_file,  mode_w, template=sam) as mul, \
             Samfile(unm_file,  mode_w, template=sam) as unm:
            for read in sam.fetch(until_eof=True):
                xa = True if len(read.tags) > 0 and read.tags[-1][0] == 'XA' else False  # have 'XA' tag or not.
                if read.mapq >= mapq_thresh and not xa:  # unique map
                    counts["unique-mapped"] += 1
                    uniq.write(read)
                elif xa:
                    counts['multiple-mapped'] += 1
                    mul.write(read)
                else:
                    counts['other'] += 1
                    unm.write(read)
    return counts


def log_split_count(pet1_count, pet2_count, log_file):
    def count_msg(counts):
        total = sum(counts.values())
        u = counts['unique-mapped']
        m = counts['multiple-mapped']
        o = counts['other']
        u_r = 0 if total == 0 else u/total
        m_r = 0 if total == 0 else m/total
        o_r = 0 if total == 0 else o/total
        msg = "unique-mapped\t{}\tpercent\t{:.2%}".format(u, u_r) + "\t"
        msg += "multiple-mapped\t{}\tpercent\t{:.2%}".format(m, m_r) + "\t"
        msg += "other\t{}\tpercent\t{:.2%}".format(o, o_r)
        return msg
    msg1 = count_msg(pet1_count)
    msg2 = count_msg(pet2_count)
    log.info("PET1 alignment results:\t" + msg1)
    log.info("PET2 alignment results:\t" + msg2)
    with open(log_file, 'w') as fo:
        fo.write("PET1\t" + msg1 + "\n")
        fo.write("PET2\t" + msg2 + "\n")


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
@click.option("--upper-tri/--non-upper-tri",
    default=True,
    help="Convert bedpe to upper-triangular format or not, default True")
@click.option("--log-file",
    default="build_bedpe.log",
    help="Sperate log file for storage some count information.")
@click.option("--bwa-log-file",
    default="bwa.log",
    help="Separate log file for storage bwa output")
def _main(file_format, input1, input2, bedpe, ncpu, bwa_index, mapq, upper_tri, log_file, bwa_log_file):
    """ Build bedpe file from fastq or sam/bam file. """
    log.info("Build bedpe from %s %s"%(input1, input2))

    outdir = dirname(bedpe)

    pre_1 = join(outdir, splitext(split(input1)[1])[0]) # prefix of input fastq file
    pre_2 = join(outdir, splitext(split(input2)[1])[0])
    if file_format == 'fastq':
        log.info("Do alignment using 'bwa aln', with {} threads".format(ncpu))
        assert bwa_index is not None
        bwa = BWA(bwa_index, algorithm='aln', log_file=bwa_log_file) # write bwa log to separate file
        bwa.run(input1, pre_1, thread=ncpu, max_diff=0, bam=True)
        bwa.run(input2, pre_2, thread=ncpu, max_diff=0, bam=True)
        bam1 = pre_1 + '.bam'
        bam2 = pre_2 + '.bam'
    else:
        # build from bam file directly
        bam1 = input1
        bam2 = input2

    log.info("split bam file.")
    with ProcessPoolExecutor(max_workers=2) as executor:
        pet1_count, pet2_count = list(executor.map(split_sam, [bam1, bam2], [mapq, mapq]))
    log_split_count(pet1_count, pet2_count, log_file)

    log.info("convert bam to bed.")
    bed1_ = pre_1 + ".bed.tmp"
    bed2_ = pre_2 + ".bed.tmp"
    p1 = subprocess.Popen("bedtools bamtobed -i {} > {}".format(pre_1+'.uniq.bam', bed1_), shell=True)
    p2 = subprocess.Popen("bedtools bamtobed -i {} > {}".format(pre_2+'.uniq.bam', bed2_), shell=True)
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