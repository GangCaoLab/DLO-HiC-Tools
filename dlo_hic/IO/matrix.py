import os
import sys
import cPickle

import numpy as np
from scipy import sparse


def save_hicmat(hicmat, file_prefix, use_sparse=True):
    """ 
    save object to npz file.
    :use_sparse: if True(default), compress array to csr sparse matrix.

    >>> save_hicmat(hicmat, "test")

    $ ls
    test.npz

    """
    prefix = file_prefix

    mat_ = hicmat.matrix

    mat = hicmat.matrix
    if use_sparse:
        mat = sparse.csr_matrix(mat)
    hicmat.matrix = None
    obj = cPickle.dumps(hicmat)
    np.savez(file_prefix + '.npz', obj=obj, matrix=mat)

    # recover matrix
    hicmat.matrix = mat_


def load_hicmat(file_name):
    """ 
    load HicMatrix/HicChrMatrix object from .npz file.

    >>> load_hicmat("test.npz")
    array([[1 2 0 ... 1]
            ...
           [1 3 0 ... 1]])
    
    """
    npz_item = np.load(file_name)
    mat = npz_item['matrix'].item()
    if type(mat) == sparse.csr.csr_matrix:
        mat = np.array(mat.todense())
    obj = npz_item['obj'].tostring()
    hicmat = cPickle.loads(obj)
    hicmat.matrix = mat
    return hicmat
