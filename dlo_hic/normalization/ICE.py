import numpy as np
from iced import filter, normalization

try:
    from dlo_hic.hicmatrix import HicChrMatrix
except ImportError:
    pass

def ice(hicmat, percentage=0.05):
    """
    ICE(iterative correction and eigenvector decomposition).
    """
    hicmat.matrix = np.array(hicmat.matrix.astype(np.float))
    if isinstance(hicmat, HicChrMatrix):
        filtered_mat = filter.filter_low_counts(hicmat.matrix,
                lengths=hicmat.lengths, percentage=percentage)
    else:
        filtered_mat = filter.filter_low_counts(hicmat.matrix,
                percentage=percentage)
    normed_mat = normalization.ICE_normalization(filtered_mat)
    hicmat.matrix = normed_mat
    hicmat.iced = True
