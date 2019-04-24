import io
import subprocess
from itertools import cycle
import logging

from dlo_hic.utils.fastaio import read_fasta, FastaRec
from dlo_hic.utils.fastqio import read_fastq
from dlo_hic.utils.suffix_tree.STree import STree


log = logging.getLogger(__name__)


N_FQ_REC = 50
SEARCH_START_POS = 76
N_BATCH = 10
BATCH_SIZE = 5
N_ALIGN = 5
MIN_LEN = 5


def multiple_alignment(input_str):
    """ Perform mafft return result fasta records """
    try:
        p = subprocess.Popen(["mafft", "-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    except FileNotFoundError:
        msg = "'mafft' for multiple alignment is not installed."
        raise FileNotFoundError(msg)
    out, err = p.communicate(input_str.encode('utf-8'))
    p.kill()
    str_io = io.StringIO(out.decode('utf-8'))
    fa_recs = [r for r in read_fasta(str_io)]
    return fa_recs


def get_fq_recs(path, n=N_FQ_REC):
    fq_recs = []
    for i, fq_rec in enumerate(read_fastq(path)):
        if i == n:
            break
        fq_recs.append(fq_rec)
    return fq_recs


def to_mafft_input(fq_recs, start_pos=SEARCH_START_POS):
    """ Convert fasta records to fasta string. for mafft input. """
    fa_recs = []
    for fq in fq_recs:
        name = fq.seqid
        seq = fq.seq
        seq = seq[start_pos-1:]
        fa = FastaRec(name, seq)
        fa_recs.append(fa)
    fa_str = "".join([str(i) for i in fa_recs])
    return fa_str


def make_batchs(fa_recs, n_batch=N_BATCH, batch_size=BATCH_SIZE):
    """ Generate multiple align result batch. """
    iter = cycle(fa_recs)
    for _ in range(n_batch):
        batch = []
        for _ in range(batch_size):
            batch.append(next(iter))
        yield batch


def find_lcs(recs):
    """ Find longest common subsequence of records. """
    st = STree([r.seq for r in recs])
    lcs = st.lcs()
    start_pos = [r.seq.find(lcs) for r in recs]
    start_pos = [p for p in start_pos if p >= 0]
    min_ = min(start_pos) if start_pos else float('inf')
    return lcs, min_


def infer_adapter_seq(fastq_path, n_fq_rec=N_FQ_REC, start_pos=SEARCH_START_POS,
                      n_batch=N_BATCH, batch_size=BATCH_SIZE, n_align=N_ALIGN, min_adapter_len=MIN_LEN):
    candidates = {}
    for _ in range(n_align):
        fq_recs = get_fq_recs(fastq_path, n_fq_rec)
        fa_str = to_mafft_input(fq_recs, start_pos)
        align_recs = multiple_alignment(fa_str)

        for recs in make_batchs(align_recs, n_batch, batch_size):
            lcs, start_pos = find_lcs(recs)
            if (len(lcs) >= min_adapter_len) and ('-' not in lcs):
                candidates.setdefault(start_pos, [])
                candidates[start_pos].append(lcs)

    if not candidates:
        return None

    candidates = list(candidates.items())
    candidates.sort(key=lambda t: len(t[1]), reverse=True)
    candidates = candidates[0][1]
    candidates.sort(key=lambda s: len(s), reverse=True)

    return candidates[0]


