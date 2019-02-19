import logging
import subprocess

import click

from dlo_hic.utils.stream import (upper_triangle, remove_redundancy,
                                  read_file, write_to_file,
                                  sort_pairs, sort_bedpe)
from dlo_hic.utils.parse_text import Bedpe, Pairs, infer_interaction_file_type


log = logging.getLogger(__name__)


@click.command(name="remove_redundancy")
@click.argument("input")
@click.argument("output")
@click.option("--distance", "-d", default=0,
    help="The threshold of distance, if pairs both ends's distance,"
         "small than this at same time, consider them as the redundancy."
         "default 0(exactly same)")
@click.option("--by-etag", "-e", default=False,
    is_flag=True,
    help="Remove redundancy by enzyme cutting site. "
         "It's useful only when input is bedpe format, and have the extends fields in file.")
@click.option("--ncpu",
    default=1,
    help="The number of cpu cores used for sort.")
def _main(input, output, distance, by_etag, ncpu):
    """
    Remove the redundancy within pairs.

    If pairs both ends's distance,
    small than the threshold distance at same time, and have same strand.
    consider them as the redundancy.

    \b
    for example:
                  +             -
        reads1    <---  ...  --->
               |-----|      |-----|
        reads2   <--- ... --->
                 +           -

        reads2 can be consider as the replection of reads1.
        reads2 will be remove.

    Arguments
    ---------
    input : str
        Path to input file, file format can be BEDPE or Pairs.

    output : str
        Path to output file.

    """

    log.info("remove redundancy on file %s"%input)

    # sort input file firstly
    log.info("transform to upper triangle form.")
    tmp0 = input + '.tmp.0'
    line_iter = read_file(input)

    fmt = infer_interaction_file_type(input)
    if fmt == Pairs:
        fmt = 'pairs'
        sort_func = sort_pairs
    else:
        fmt = 'bedpe'
        sort_func = sort_bedpe

    line_iter = upper_triangle(line_iter, fmt=fmt)
    write_to_file(line_iter, tmp0)

    log.info("sorting and remove redundancy ...")
    line_iter = sort_func(tmp0, ncpu=ncpu)
    line_iter = remove_redundancy(line_iter, fmt, distance, by_etag=by_etag)
    write_to_file(line_iter, output)

    subprocess.check_call(['rm', tmp0])  # remove tmp files
    log.info("Remove redundancy done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")