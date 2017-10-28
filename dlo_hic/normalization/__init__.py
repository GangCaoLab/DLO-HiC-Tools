import numpy as np

from .ICE import ice


def remove_zero(matrix, method='small'):
    """
    remove zero values in hicmat.matrix.
    :method: 'min', replace zero by min val.
             'max', replace zero by max value.
             'tiny': replace zero by smallest float value.
             'small':  replace zero by a very small float value.(default)
    """
    import dlo_hic.hicmatrix
    tiny = np.finfo(np.float).tiny
    small = np.float(1e-2)

    if isinstance(matrix, dlo_hic.hicmatrix.HicMatrix):
        mat = matrix.matrix
    else:
        mat = matrix

    if method == 'min':
        min_val = mat[hicmat.matrix > 0].min()
        mat[mat == 0] = min_val
    if method == 'max':
        max_val = mat.max()
        mat[mat == 0] = max_val
    if method == 'tiny':
        mat[mat == 0] = tiny
    if method == 'small':
        mat[mat == 0] = small
