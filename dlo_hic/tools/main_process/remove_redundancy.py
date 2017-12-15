import sys
import subprocess

import click

from dlo_hic.utils import read_args
from dlo_hic.utils.wrap.tabix import sort_bedpe_reads1
from dlo_hic.utils.parse_text import Bedpe


def bedpe_upper_triangle(bedpe_file, output):
    """
    transform bedpe file's all line to upper trangle form.
    """
    with open(bedpe_file) as fi, open(output, 'w') as fo:
        for line in fi:
            bpe = Bedpe(line)
            bpe.to_upper_trangle()
            outline = str(bpe) + "\n"
            fo.write(outline)


def remove_redundancy(bedpe_file, output, distance):
    """
    remove redundancy lines in the bedpe file,
    input file must in upper trangle form and sorted according to reads1(first 3 cols).
    """
    with open(bedpe_file, 'r') as f, open(output, 'w') as fo:
       base = Bedpe(f.readline())
       while True:
            for line in f:
                another = Bedpe(line)
                if base.is_rep_with(another, distance): # is replication, check next line.
                    continue
                else: # not replication, output base line and change base line.
                    out_line = str(base) + "\n"
                    fo.write(out_line)
                    base = another 
                    break
            else: # arrive at end of file.
                out_line = str(base) + "\n"
                fo.write(out_line)
                break


@click.command(name="remove_redundancy")
@click.argument("bedpe")
@click.argument("output")
@click.option("--distance", "-d", default=50,
    help="The threshold of distance, if pairs both ends's distance,"
         "small than this at same time, consider them as the redundancy.")
def _main(bedpe, output, distance):
    """
    Remove the redundancy within pairs.
    
    If pairs both ends's distance,
    small than the threshold distance at same time,
    consider them as the redundancy.

    \b
    for example:
        reads1    <---  ...  --->
               |-----|      |-----|
        reads2   <--- ... --->

        reads2 can be consider as the replection of reads1.
        reads2 will be remove.


    """

    # sort input file firstly
    tmp0 = bedpe + '.tmp.0'
    bedpe_upper_triangle(bedpe, tmp0)
    tmp1 = bedpe + '.tmp.1'
    sort_bedpe_reads1(tmp0, tmp1)
    remove_redundancy(tmp1, output, distance)

    subprocess.check_call(['rm', tmp0]) # remove tmp files
    subprocess.check_call(['rm', tmp1])


main = _main.callback


if __name__ == "__main__":
    _main()