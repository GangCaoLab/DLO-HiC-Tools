import io
import subprocess
from itertools import cycle
import logging
from collections import OrderedDict, defaultdict

from dlo_hic.utils.fastaio import read_fasta, FastaRec
from dlo_hic.utils.fastqio import read_fastq
from dlo_hic.utils.suffix_tree.STree import STree


log = logging.getLogger(__name__)


N_FQ_REC = 50
SEARCH_START_POS = 0
PROB_THRESH = 0.75


def multiple_alignment(input_str, result_file=None):
    """ Perform `mafft` return result fasta records """
    try:
        p = subprocess.Popen(["mafft", "-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    except FileNotFoundError:
        msg = "'mafft' for multiple alignment is not installed."
        raise FileNotFoundError(msg)
    out, err = p.communicate(input_str.encode('utf-8'))
    p.kill()
    str_io = io.StringIO(out.decode('utf-8'))
    fa_recs = [r for r in read_fasta(str_io)]
    if result_file is not None:
        with open(result_file, 'w') as f:
            f.write(out.decode('utf-8'))
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
        seq = seq[start_pos:]
        fa = FastaRec(name, seq)
        fa_recs.append(fa)
    fa_str = "".join([str(i) for i in fa_recs])
    return fa_str


def find_lcs(recs):
    """ Find longest common subsequence of records. """
    st = STree([r.seq for r in recs])
    lcs = st.lcs()
    start_pos = [r.seq.find(lcs) for r in recs]
    start_pos = [p for p in start_pos if p >= 0]
    min_ = min(start_pos) if start_pos else float('inf')
    return lcs, min_


def fa_to_seqs(fa_recs):
    return [rec.seq for rec in fa_recs]


def get_char_prob(seqs):
    """ count char(base) occurace probability at each position. """
    seq_len = len(seqs[0])
    probs = []
    for i in range(seq_len):
        p = {c:0 for c in 'ATCG-'}
        for s in seqs:
            c = s[i].upper()
            if c != 'N':
                p[c] += 1
        probs.append(p)
    probs = [{k:v/sum(p.values()) for k,v in p.items()} for p in probs]
    return probs


def plot_char_prob_stacked_bar(probs, start_pos=0):
    """ stacked bar plot of char probability """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import pandas as pd

    base2color = OrderedDict({
        'A': '#b5ffb9', 
        'T': '#f9bc86',
        'C': '#a3acff',
        'G': '#ff9c9c',
        '-': 'grey',
    })
    
    r = [start_pos + i for i in range(len(probs))]

    raw_data = defaultdict(list)
    for p in probs:
        for b, v in p.items():
            raw_data[b].append(v)
    df = pd.DataFrame(raw_data)
    
    # plot
    fig, ax = plt.subplots()
    # Create Bars
    bottom = 0
    for b, color in base2color.items():
        ax.bar(r, df[b], bottom=bottom, color=color, edgecolor='white')
        bottom = (bottom + df[b]) if (bottom is not 0) else df[b]
    # modify
    plt.xlabel("Positions in mafft alignment result")
    plt.ylabel("Percentage")
    plt.legend(handles=[Patch(color=c, label=l) for l,c in base2color.items()])
    plt.xlim(start_pos-1, start_pos+len(probs))
    plt.ylim(0, 1)
    return fig, ax



def infer_conserved(probs, thresh=PROB_THRESH):
    chars = []
    for p in probs:
        for b, v in p.items():
            if v >= thresh:
                chars.append(b)
                break
        else:
            chars.append('N')
    return "".join(chars)


def infer_adapter(fq_path, n_fq_rec=N_FQ_REC, search_start_pos=SEARCH_START_POS, prob_thresh=PROB_THRESH):
    """ inference adapter sequence from a fastq file
    
    Arguments
    ---------
    fq_path : str
        Path to input Fastq file.

    n_fq_rec : int
        Number of fastq records used for inference adapter sequence.

    search_start_pos : int
        Start position of adapter search.

    prob_thresh : float
        Threshold of occurace probability in adapter inference.

    Return
    ------
    probs : [dict]
        char(base) occurace probability at each position.

    flag_seq : str
        flag sequence used for inference adapter sequence.
    
    adapter_seq : str
        adapter sequence.
    """
    mafft_in = to_mafft_input(get_fq_recs(fq_path, n_fq_rec), search_start_pos)
    fa_recs = multiple_alignment(mafft_in)
    probs = get_char_prob(fa_to_seqs(fa_recs))
    flag_seq = infer_conserved(probs, prob_thresh)
    adapter_seq = sorted([s.replace('-', '') for s in flag_seq.split("N")], key=lambda s: len(s))[-1]
    return probs, flag_seq, adapter_seq


def infer_adapter_seq(fq_path, n_fq_rec=N_FQ_REC, search_start_pos=SEARCH_START_POS, prob_thresh=0.5):
    return infer_adapter(fq_path, n_fq_rec, search_start_pos, prob_thresh)[2]

