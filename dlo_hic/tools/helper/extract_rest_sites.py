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
import h5py
import numpy as np

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

        out_chunk = [] 
        for match in re.finditer(rest, seq, re.IGNORECASE):
            out_chunk.append(
                match.start()
            )
        output_queue.put( (chr_, '+', out_chunk) )

        if rest_rc != rest: # find reverse complement restriction site
            out_chunk = []
            for match in re.finditer(rest_rc, seq, re.IGNORECASE):
                out_chunk.append(
                    match.start()
                )
            output_queue.put( (chr_, '-', out_chunk) )


def outputer(output, out_fmt, rest_seq, output_queue):
    """ output extracted results """
    if out_fmt == 'tab':
        with open(output, 'w') as f:
            while 1:
                out_tupl = output_queue.get()
                if out_tupl is None:
                    log.debug("Process-output done.")
                    f.flush()
                    break
                chr_, strand, out_chunk = out_tupl
                for start in out_chunk:
                    # BED6: [chr, start, end, name, score, strand]
                    out_itms = [chr_, str(start), str(start + len(rest_seq)), ".", "0", strand]
                    out_line = "\t".join(out_itms) + "\n"
                    f.write(out_line)
    elif out_fmt == 'hdf5':
        output = output + '.hdf5' if not output.endswith('.hdf5') else output
        with h5py.File(output, 'w') as f:
            f.create_group("chromosomes")
            f.attrs['rest_seq'] = rest_seq
            while 1:
                out_tupl = output_queue.get()
                if out_tupl is None:
                    log.debug("Process-output done.")
                    f.flush()
                    break
                chr_, _, out_chunk = out_tupl
                f.create_dataset("chromosomes/"+chr_, data=np.array(out_chunk, dtype=np.int))
    else:
        raise NotImplementedError("output format only support tab and hdf5.")


@click.command(name="extract_rest_sites")
@click.argument("fasta", nargs=1)
@click.option("--rest-seq", "-r", "rest", required=True,
    help="The sequence of restriction site")
@click.argument("output", nargs=1)
@click.option("--output-format", "-f",
    type=click.Choice(['tab', 'hdf5']),
    default='tab',
    show_default=True,
    help="The file output format.")
@click.option("--processes", "-p", default=1,
    help="Use how many processes to run. default 1")
def _main(fasta, rest, output, output_format, processes):
    """ Extract all restriction sites from fasta file, save to BED6 file format. """

    log.info("Extract restriction sites %s from %s"%(rest, fasta))

    if output.endswith(".gz"):
        output = output.replace(".gz", "")

    faidx = pyfaidx.Fasta(fasta)
    chrs = faidx.keys()

    task_queue   = mp.Queue()
    output_queue = mp.Queue(maxsize=processes)

    processes = min(processes, len(chrs))

    workers = [mp.Process(target=worker,
                       args=(task_queue, output_queue, rest, fasta))
               for i in range(processes)]

    for chr_ in chrs:
        task_queue.put(chr_)
    for w in workers: # put end flag to task queue
        task_queue.put(None)

    output_p = mp.Process(target=outputer, args=(output, output_format, rest, output_queue))

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

    if output_format == 'tab':
        # sort output bed file
        log.info("sorting bed file ...")
        sort_bed6(output, output+'.sort')
        subprocess.check_call(["mv", output+'.sort', output])
        log.info("building tabidx...")
        index_bed6(output)
        subprocess.check_call(["rm", output])
        log.info("Result storaged in bgziped file %s"%(output+".gz"))

main = _main.callback

if __name__ == "__main__":
    eval("_main()")