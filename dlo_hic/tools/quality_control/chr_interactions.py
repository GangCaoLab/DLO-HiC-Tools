import logging
import subprocess as subp

import numpy as np
import pandas as pd
import click

log = logging.getLogger(__name__)

def get_counter(path_pairs):
    def counter(chr1, chr2):
        cmd = "pairix {} '{}|{}' | wc -l".format(path_pairs, chr1, chr2)
        log.debug(cmd)
        logging
        p = subp.Popen(cmd, stdout=subp.PIPE, shell=True)
        out, _ = p.communicate()
        num = int(out.decode('utf-8'))
        return num
    return counter

def get_chromosomes(chromosome_file):
    chrs = []
    with open(chromosome_file) as f:
        for line in f:
            chrs.append(line.strip().split()[0])
    return chrs

def fetch_chr_pairs(path_pairs):
    cmd = "pairix -l {}".format(path_pairs)
    p = subp.Popen(cmd, shell=True, stdout=subp.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        chr1, chr2 = line.strip().split('|')
        yield (chr1, chr2)

def chromosome_combinations(path_pairs, chromosome_file):
    possible_chrs = set(get_chromosomes(chromosome_file))
    for chr1, chr2 in fetch_chr_pairs(path_pairs):
        if (chr1 not in possible_chrs) or (chr2 not in possible_chrs):
            continue
        yield (chr1, chr2)

def count_all(path_pairs, chromosome_file):
    chrs = get_chromosomes(chromosome_file)
    chrs_num = len(chrs)
    counter = get_counter(path_pairs)
    mat = np.zeros((chrs_num, chrs_num), dtype=np.int)
    df = pd.DataFrame(mat, columns=chrs, index=chrs)
    for chr1, chr2 in chromosome_combinations(path_pairs, chromosome_file):
        count = counter(chr1, chr2)
        df.loc[chr1, chr2] = count
        df.loc[chr2, chr1] = count
    return df

def save_to_file(df, output):
    df.to_csv(output)


@click.command(name="chr_interactions")
@click.argument("input")
@click.argument("chromosome-file")
@click.argument("output")
def _main(input, chromosome_file, output):
    """
    Count the interactions between chromosomes.

    \b
    Input:
        1. Indexed pairs gz(bgzip) file
        2. Chromosome file:
            A text file contain chromosome's name at the first column.
            like:
            \b
            \"\"\"
            chr1
            chr2
            chr3
            chr4
            chr5
            ...
            \"\"\"
    
    \b
    Output:
        A csv file contain a interaction counts matrix between chromosomes.

    """
    df = count_all(input, chromosome_file)
    save_to_file(df, output)

main = _main.callback


if __name__ == "__main__":
    eval("_main()")
