"""
Command line tool for build bedpe file from fastq or sam/bam file. 
"""

import argparse

import os
import subprocess
import pybedtools
import multiprocessing

from dlo_hic.utils import BWA
from dlo_hic.utils import read_args

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("input1")
    parser.add_argument("input2")
    parser.add_argument("output")

    parser.add_argument("--file-format", "-f", dest="file_format", required=True,
        choices=['fastq', 'bam', 'sam'],
        help="The file type of input file, if 'fastq' will preform bwa alignment firstly.")

    parser.add_argument("--threads", "-t",
        type=int,
        default=multiprocessing.cpu_count(),
        help="how mant threads used to run bwa")

    parser.add_argument("--bwa-index", dest="bwa_index",
        help="The bwa index prefix.")

    parser.add_argument("--mapq",
        type=int,
        default=20,
        help="the mapq threshold used to filter mapped records.")

    return parser


def beds2bedpe(bed1, bed2, bedpe_filename):
    # build bed's index
    name2interval = {}
    for interval in bed1[:]:
        if interval.name in name2interval:
            raise ValueError("one reads occurred twice")
        name2interval[interval.name] = interval
    with open(bedpe_filename, 'w') as f:
        for i2 in bed2[:]:
            i1 = name2interval.get(i2.name, None)
            if i1:
                fields = [i1.chrom, i1.start, i1.end, i2.chrom, i2.start, i2.end,
                          i1.name, 1, i1.strand, i2.strand]
                fields = [str(i) for i in fields]
                line = "\t".join(fields) + "\n"
                f.write(line)


def main(file_format, input1, input2, output, threads, bwa_index, mapq):
    pre_1 = os.path.splitext(input1)[0] # prefix
    pre_2 = os.path.splitext(input2)[0]
    if file_format == 'fastq':
        # alignment firstly if the input format is fastq
        assert bwa_index is not None
        bwa = BWA(bwa_index)
        bwa.run(input1, pre_1, thread=threads, mem=False)
        bwa.run(input2, pre_2, thread=threads, mem=False)
        sam1 = pre_1 + '.sam'
        sam2 = pre_2 + '.sam'
    else:
        sam1 = input1
        sam2 = input2
    subprocess.check_call("samtools view {} -b -q {} > {}.filtered.bam".format(sam1, mapq, pre_1), shell=True)
    subprocess.check_call("samtools view {} -b -q {} > {}.filtered.bam".format(sam2, mapq, pre_2), shell=True)
    bed1 = pybedtools.BedTool(pre_1+".filtered.bam").bam_to_bed().saveas(pre_1+".bed")
    bed2 = pybedtools.BedTool(pre_2+".filtered.bam").bam_to_bed().saveas(pre_2+".bed")
    beds2bedpe(bed1, bed2, output)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()
    read_args(args, globals())
    main(file_format, input1, input2, output, threads, bwa_index, mapq)