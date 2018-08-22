import io
import gzip
import logging
from collections import namedtuple

import click

from dlo_hic.utils.parse_text import Bedpe, Pairs, Bedpe_err


log = logging.getLogger(__name__)


def log_counts(counts, log_file=None):

    def cal_ratio(c, total):
        try:
            return c / total
        except ZeroDivisionError:
            return 0

    total = counts.total
    inter = counts.inter
    intra = counts.intra
    long_range = counts.long_range
    if not log_file:
        log.info("Interaction Quality Control:")
        log.info("total: %d"%total)

        log.info("\t%s\t%d\t%.5f"%('inter-chromosome', inter, cal_ratio(inter, total)))
        log.info("\t%s\t%d\t%.5f"%('intra-chromosome', intra, cal_ratio(intra, total)))
        log.info("\t%s\t%d\t%.5f"%('long-range', long_range, cal_ratio(long_range, intra)))

    else:
        with open(log_file, 'w') as f:
            f.write("total\t{}\n".format(total))
            f.write("%s\t%d\n"%('inter-chromosome', inter))
            f.write("%s\t%d\n"%('intra-chromosome', intra))
            f.write("%s\t%d\n"%('long-range', long_range))


class InteractionCounts(object):
    def __init__(self):
        self.inter = 0
        self.intra = 0
        self.long_range = 0
        self.total = 0
        self.abnormal_1 = 0
        self.abnormal_2 = 0


def open_file(file_name):
    if file_name.endswith(".gz"):
        f = io.TextIOWrapper(gzip.open(file_name))
    else:
        f = open(file_name)
    return f


@click.command(name="interactions_qc")
@click.argument("input")
@click.option("--log-file",
    default="interactions_qc.txt",
    help="Sperate to store qc information. default: interactions_qc.txt")
@click.option("--long-range-cutoff",
    default=20*(10**3),
    help="The cutoff of intra long range interaction, "
         "span larger than this count as long-range. default 20Kb")
def _main(input, log_file, long_range_cutoff):
    """
    Count ratio of:

    \b
        inter-chromosome
        intra-chromosome
        long-range

    in interactions data (in bedpe or pairs format, support gziped file).
    """

    # init counter dict
    counter = InteractionCounts()

    # conform input file format
    if input.endswith("pairs") or input.endswith("pairs.gz"):
        fmt = Pairs
        log.info("input pairs file.")
    elif input.endswith("bedpe") or input.endswith("bedpe.gz"):
        fmt = Bedpe
        log.info("input bedpe file.")
    elif input.endswith("bedpe.err") or input.endswith("bedpe.err.gz"):
        # count abnormal bedpe file(produced by noise reduce)
        fmt = Bedpe_err
        log.info("input abnormal bedpe file")
    else:
        raise NotImplementedError("Only support pairs and bedpe file format.")

    with open_file(input) as f:
        for line in f:
            pair = fmt(line)
            if pair.chr1 != pair.chr2:
                counter.inter += 1
            else:
                counter.intra += 1
                span = abs(pair.pos2 - pair.pos1)
                if span >= long_range_cutoff:
                    counter.long_range += 1

            counter.total += 1

            if fmt is Bedpe_err:
                if pair.abormal_type == 'abnormal-1':
                    counter.abnormal_1 += 1
                else:
                    counter.abnormal_2 += 1

    log_counts(counter)
    log_counts(counter, log_file=log_file)


main = _main.callback


if __name__ == "__main__":
    _main()
