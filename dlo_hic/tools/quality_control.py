# -*- coding: utf-8 -*-

"""
quality_control.py
~~~~~~~~~~~~~~~~~~

Generate the quality control report of dlo-hic pipeline.

quality control indexes of one:
+---+--------------------------------------+----------------------+
|id | BASIC ITEMS                          | file partten         |
+---+--------------------------------------+----------------------+
| 1 | raw reads                            |                      |
+---+--------------------------------------+----------------------+
| 2 | clean reads                          | clean_.*\.fq         |
+---+--------------------------------------+----------------------+
| 3 | linker reads                         | hookers_.*\.fq       |
|   |     A-B B-A linker                   |                      |
|   |     B-B A-B linker (optional, TODO)  |                      |
+---+--------------------------------------+----------------------+
| 4 | mapping reads                        | mapping_.*\.sam      |
|   |     unique mapping reads             | filtered_.*\.bam     |
|   |     paired unique mapping read       | pair_.*\.bedpe       |
+--++--------------------------------------+----------------------+
| 5 | non-redundant mapped reads           | nr_pair_.*\.bedpe    |
+---+--------------------------------------+----------------------+
| 6 | final reads (after noise reduction.) | clean_pair_.*\.bedpe |
+---+--------------------------------------+----------------------+

"""

import re
import os
import sys
import json
import argparse
from subprocess import Popen, PIPE

import pdb


def argument_parser():
    parser = argparse.ArgumentParser(
            description = "Creat interaction matrix from pairs file in bedpe format.")

    parser.add_argument("--dir", "-d",
            required=True,
            help="The directory of the dlo-hic pipeline intermediate files.")

    parser.add_argument("--raw", "-r",
            required=True,
            help="raw data, the input of pipeline.")

    parser.add_argument("--output", "-o",
            type=argparse.FileType(mode="w"),
            default=sys.stdout,
            help="output file, default stdout.")

    parser.add_argument("--format", "-f",
            default="tsv",
            help="quality control report file format,"
                 "options: json, tsv, csv"
                 "default in json format.")

    return parser


def items():
    """ generate a dict store name parttens of target files. """
    item2partten = {
            0:   r".*\.(fq|sam|bam|bedpe)$", # for grep all considered files
            2.0: r"clean.*\.fq$",
            3.0: r"hookers_.*\.fq$",
            4.0: r"mapping_.*\.sam$",
            4.1: r"filtered_.*\.bam$",
            4.2: r"pair_.*\.bedpe$",
            5.0: r"nr_pair_.*\.bedpe$",
            6.0: r"clean_pair_.*\.bedpe$",
            }
    item2name = {
            1.0: "raw reads",
            2.0: "clran reads",
            3.0: "linker reads",
            4.0: "mapped reads",
            4.1: "unique mapped reads",
            4.2: "paired unique mapped reads",
            5.0: "non redundant reads",
            6.0: "noise reduced reads"
            }
    return item2partten, item2name


def count_reads(files):
    """ count all file's reads number """

    def gen_counter(prefix=None):
        """ high order function, generate a counter function """
        def counter_process(fname):
            if prefix:
                return Popen("cat {fname} | {command} | wc -l ".format(command=prefix, fname=fname),
                        shell=True, stdout=PIPE)
            else:
                return Popen("wc -l {}".format(fname), shell=True, stdout=PIPE)
        return counter_process

    text_counter_process = gen_counter()
    sam_counter_process = gen_counter("grep -v '^@'")
    bam_counter_process = gen_counter("samtools view")

    # generate counters
    counters = []
    for fname in files:
        if fname.endswith(".sam"):
            counters.append(sam_counter_process(fname))
        elif fname.endswith(".bam"):
            counters.append(bam_counter_process(fname))
        else:
            counters.append(text_counter_process(fname))

    fetch_line_num = lambda p: int(p.communicate()[0].split()[0])

    # count all files
    counts = []
    for fname, counter in zip(files, counters):
        if fname.endswith(".fq") or fname.endswith(".fastq"):
            counts.append(fetch_line_num(counter)/4) # fastq file 4 lines as one reads
        else:
            counts.append(fetch_line_num(counter))

    return counts


def cal_ratios(counts):
    """ calculate count ratios """
    max_ = max(counts)
    max_ = float(max_)
    return [count/max_ for count in counts]


def gen_report(infos, format):
    """ generate qc report. """
    infos = sorted(infos, key=lambda pair:pair[0])

    def gen_sv(s):
        report = ""
        for items in infos:
            items = [str(i) for i in items]
            report += (s.join(items) + "\n")
        return report

    if format == "tsv":
        report =  gen_sv(s="\t")
    elif format == "csv":
        report =  gen_sv(s=",")
    else:
        raise NotImplementedError("report format is not support.")

    return report


def main(items=items):
    parser = argument_parser()
    args = parser.parse_args()
    assert os.path.isdir(args.dir)
    assert os.path.isfile(args.raw)
    assert args.raw.endswith(".fq") or args.raw.endswith(".fastq")
    item2partten, item2name = items()
    partten_all = item2partten.pop(0)
    files = [os.path.join(args.dir, f) # fetch all considered files
            for f in os.listdir(args.dir) if re.match(partten_all, f)]
    files.append(args.raw)
    # classify files
    def classify(fname):
        for i, partten in item2partten.items():
            if re.match(partten, os.path.basename(fname)):
                return i
        return 1.0
    class_ = [classify(fname) for fname in files]
    names = [item2name[i].replace(" ", "_") for i in class_]
    reads_nums = count_reads(files) # count reads of all files
    ratios = cal_ratios(reads_nums)
    infos = zip(class_, names, files, reads_nums, ratios)
    report = gen_report(infos, args.format) # generate qc report
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
