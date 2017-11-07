from copy import copy

import numpy as np

import dlo_hic.hicmatrix

def call_diff(hicmatA, hicmatB):
    """ call difference(log2 FoldChange) between hicmatA and hicmatB. """

    matA = np.copy(hicmatA.matrix)
    matB = np.copy(hicmatB.matrix)
    assert matA.shape == matB.shape

    # remove zero values in matrix
    min_val = min(matA[matA > 0].min(), matB[matB > 0].min()) 
    matA[matA == 0] = min_val
    matB[matB == 0] = min_val

    diff = matB / matA # fold change
    diff = np.log2(diff) # log2 FC
    diff = diff - diff.mean()
    diff[diff == np.inf] = diff[diff != np.inf].max() # remove inf value
    diff[diff == -np.inf] = diff[diff != -np.inf].min() # remove -inf value
    diff[np.isnan(diff)] = 0 # remove nan value

    diff_mat = copy(hicmatA)
    diff_mat.matrix = diff
    diff_mat.__class__ = dlo_hic.hicmatrix.DiffMatrix
    return diff_mat