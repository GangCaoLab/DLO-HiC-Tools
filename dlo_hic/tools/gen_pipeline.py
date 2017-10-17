# -*- coding: utf-8 -*-

"""
gen_pipeline.py
~~~~~~~~~~~~~~~
Generate a DLO-Hic analyze pipeline(bash script).

"""

import os
import sys
import argparse
import subprocess


file_dir = os.path.dirname(os.path.realpath(__file__))
template_dir = os.path.join(file_dir, "../template")

default_temp    = os.path.join(template_dir, "pipeline.txt")
example_linkers = os.path.join(template_dir, "linkers.json")

example_chromosome_length_GRCG38 = os.path.join(template_dir, "example_chromosome_length_GRCG38")
example_chromosome_length_GRCG37 = os.path.join(template_dir, "example_chromosome_length_GRCG37")
example_chromosome_length_hg38 = os.path.join(template_dir, "example_chromosome_length_hg38")
example_chromosome_length_hg19 = os.path.join(template_dir, "example_chromosome_length_hg19")


def argument_parser():
    parser = argparse.ArgumentParser(description=
            "Generate a DLO-Hic analyze pipeline(bash script).")
    
    parser.add_argument("--template",
            default=default_temp,
            help="The pipeline template file.")

    parser.add_argument("--example_linkers",
            help="generate a example of linkers config json file.")

    parser.add_argument("--referece", "-r",
            default="GRCH38",
            help="The reference used, option: GRCH38, GRCH37, hg38, hg19"
                 "default: GRCH38")

    parser.add_argument("-Q",
            default="64",
            help="The quality standard of fastq. default 64")

    parser.add_argument("--preprocess", "-p",
            action="store_true",
            help="do preprocess(remove P7 adapter) if specified.")

    parser.add_argument("--adapter",
            help="P7 adapter sequence.")

    parser.add_argument("--trim_barcode", "-t",
            action="store_true",
            help="trim barcode if specified.")

    parser.add_argument("--barcode_left", "-bl",
            type=int,
            help="trim how many length base at left side.")

    parser.add_argument("--barcode_right", "-br",
            type=int,
            help="trim how many length base at right side.")

    parser.add_argument("--thread_num", "-@",
            type=int,
            help="using how many cpu cores.")

    parser.add_argument("--sub_part_size",
            type=int,
            default=10000,
            help="sub part size in parallel process. This is not very important,"
            " maybe you just use default value, that is ok, default 10000")

    return parser


def substitution_table(qua, prep, adapter=None):
    """ generate a substitution_table for render template. """
    sub_table = {
            "BATH_PATH": subprocess.check_call("which bash", shell=True).strip(),
            "QUA":               qua,
            "PREP":              prep,
            "ADAPTER":           adapter,
            "TRIM_BARCODE":      trim_barcode,
            "BARCODE_LEFT":      barcode_left,
            "BARCODE_RIGHT":     barcode_right,
            "LINKERS_CONFIG":    linker_config,
            "CHROMOSOME_LENGTH": chromosome_length,
            "THREAD_NUM":        thread_num,
            "SUB_PART_SIZE":     sub_part_size,
            }
    return sub_table


def main():
    parser = argument_parser()
    args = parser.parse_args()


if __name__ == "__main__":
    main()
