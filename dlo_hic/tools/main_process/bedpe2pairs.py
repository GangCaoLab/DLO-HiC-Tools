import logging
import subprocess

import click

from dlo_hic.utils.wrap.tabix import sort_pairs, index_pairs
from dlo_hic.utils.parse_text import Bedpe

log = logging.getLogger(__name__)


def bedpe2pairs(input, output):
    with open(input) as fi, open(output, 'w') as fo:
        for line in fi:
            bpe = Bedpe(line)
            pairs_line = bpe.to_pairs_line()
            fo.write(pairs_line + "\n")


@click.command(name="bedpe2pairs")
@click.argument("bedpe", nargs=1)
@click.argument("pairs", nargs=1)
@click.option("--keep/--no-keep", default=False,
    help="keep non compressed pairs file, " + \
         "if you need create .hic file use --keep option.")
def _main(bedpe, pairs, keep):
    """
    Transform bedpe format file to pairs format, and index it use pairix

    \b
    about pairs format:
    https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md
    """
    log.info("convert %s to pairs file %s ..."%(bedpe, pairs))
    # sort input file firstly
    tmp0 = bedpe + '.tmp.0'
    bedpe2pairs(bedpe, tmp0)
    log.info("sorting pairs ...")
    # add header
    header = "## pairs format v1.0\n" +\
             "#columns: readID chr1 position1 chr2 position2 strand1 strand2\n"
    with open(pairs, 'w') as f:
        f.write(header)
    sort_pairs(tmp0, pairs)
    subprocess.check_call(['rm', tmp0]) # remove tmp files
    index_pairs(pairs)
    if keep:
        log.info("pairs file with header storaged at %s"%pairs)
        log.info("pairix indexed bgziped pairs file storaged at %s"%(pairs+'.gz'))
    else:
        log.info("pairix indexed bgziped pairs file storaged at %s"%(pairs+'.gz'))
        subprocess.check_call(['rm', pairs]) # remove uncompressed file


main = _main.callback


if __name__ == "__main__":
    _main()