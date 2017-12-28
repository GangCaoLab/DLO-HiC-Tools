import io
import gzip
import logging

import click

from dlo_hic.utils.parse_text import Bedpe, Pairs


log = logging.getLogger(__name__)


def log_counts(counts, log_file=None):
    if not log_file:
        log.info("Interaction Quality Control:")
        total = counts.pop('total')
        log.info("total: %d"%total)
        for k, v in counts.items():
            try:
                ratio = v / total
            except ZeroDivisionError:
                ratio = 0
            log.info("\t%s\t%d\t%.5f"%(k, v, ratio))
    else:
        with open(log_file, 'w') as f:
            for k, v in counts.items():
                outline = "\t".join([str(k), str(v)]) + "\n"
                f.write(outline)


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
@click.option("--file-format",
    default="pairs",
    type=click.Choice(["pairs", "bedpe"]),
    help="The file format of interactions data. default pairs format")
@click.option("--long-range-cutoff",
    default=20*(10**3),
    help="The cutoff of intra long range interaction, "
         "span larger than this count as long-range. default 20Kb")
def _main(input, log_file, file_format, long_range_cutoff):
    """
    Count ratio of:

    \b
        inter-chromosome
        intra-chromosome
        long-range

    in interactions data (in bedpe or pairs format, support gziped file).
    """
    # init counter dict
    counter = {
        "total": 0,
        "inter-chromosome": 0,
        "intra-chromosome": 0,
        "long-range": 0,
    }

    fmt = Pairs if file_format == 'pairs' else Bedpe

    with open_file(input) as f:
        for line in f:
            pair = fmt(line)
            if pair.chr1 != pair.chr2:
                counter["inter-chromosome"] += 1
            else:
                counter["intra-chromosome"] += 1
                span = abs(pair.pos2 - pair.pos1)
                if span >= long_range_cutoff:
                    counter["long-range"] += 1

            counter["total"] += 1

    log_counts(counter)
    log_counts(counter, log_file=log_file)


main = _main.callback


if __name__ == "__main__":
    _main()
