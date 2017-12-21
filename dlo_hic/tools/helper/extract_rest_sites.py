import re
import sys
import time
import tempfile
import logging
import subprocess
from queue import Empty
import multiprocessing as mp

import pyfaidx
import click

from dlo_hic.utils import read_args
from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils.wrap.tabix import sort_bed6, index_bed6


log = logging.getLogger(__name__)


def worker(task_queue, output_queue, rest, fasta):
    rest = str(rest)
    rest_rc = rc(rest)
    faidx = pyfaidx.Fasta(fasta)
    while 1:
        chr_ = task_queue.get()
        if chr_ is None:
            log.debug("Process-%d done"%mp.current_process().pid)
            break

        seq = faidx[chr_][:].seq # read sequence

        for match in re.finditer(rest, seq, re.IGNORECASE):
            output_queue.put(
                (chr_, str(match.start()), str(match.end()), '.', '0', '+')
            )

        if rest_rc != rest: # find reverse complement restriction site
            for match in re.finditer(rest_rc, seq, re.IGNORECASE):
                output_queue.put(
                    (chr_, str(match.start()), str(match.end()), '.', '0', '-')
                )


def outputer(output_file, output_queue):
    """ output extracted results """
    while 1:
        record = output_queue.get()
        if record is None:
            log.debug("Process-output done.")
            output_file.flush()
            break
        line = "\t".join(record) + "\n"
        output_file.write(line)


@click.command(name="extract_rest_sites")
@click.argument("fasta", nargs=1)
@click.option("--rest-seq", "-r", "rest", required=True,
    help="The sequence of restriction site")
@click.argument("output", nargs=1)
@click.option("--processes", "-p", default=1,
    help="Use how many processes to run. default 1")
def _main(fasta, rest, output, processes):
    """ Extract all restriction sites from fasta file, save to BED6 file format. """

    log.info("Extract restriction sites %s from %s"%(rest, fasta))

    if output.endswith(".gz"):
        output = output.replace(".gz", "")

    faidx = pyfaidx.Fasta(fasta)
    chrs = faidx.keys()

    task_queue   = mp.Queue()
    output_queue = mp.Queue()

    processes = min(processes, len(chrs))

    workers = [mp.Process(target=worker,
                       args=(task_queue, output_queue, rest, fasta))
               for i in range(processes)]

    for chr_ in chrs:
        task_queue.put(chr_)
    for w in workers: # put end flag to task queue
        task_queue.put(None)

    with tempfile.NamedTemporaryFile(mode='w') as tmp:
        output_p = mp.Process(target=outputer, args=(tmp, output_queue))

        # launch processes
        for w in workers:
            w.start()
        output_p.start()

        log.info("%d workers spawned for extract restriction sites"%len(workers))

        # wait all tasks done
        for w in workers:
            w.join()

        # put end flag to output queue
        output_queue.put(None)

        output_p.join() # wait all result outputed

        # sort output bed file
        log.info("sorting bed file ...")
        sort_bed6(tmp.name, output)
        log.info("building tabidx...")
        index_bed6(output)
        subprocess.check_call(["rm", output])
        log.info("Result storaged in bgziped file %s"%(output+".gz"))

main = _main.callback

if __name__ == "__main__":
    _main()