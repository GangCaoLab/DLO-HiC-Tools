import numpy as np

from .ICE import ice

tiny = np.finfo(np.float).tiny
small = np.float(1e-10)

def remove_zero(hicmat, method='tiny'):
    """
    remove zero values in hicmat.matrix.
    :method: 'min', replace zero by min val.
             'max', replace zero by max value.
             'tiny': replace zero by smallest float value.
             'small':  replace zero by a very small float value.(default)
    """
    if method == 'min':
        min_val = hicmat.matrix[hicmat.matrix > 0].min()
        hicmat.matrix[hicmat.matrix == 0] = min_val
    if method == 'max':
        max_val = self.matrix.max()
        hicmat.matrix[hicmat.matrix == 0] = max_val
    if method == 'tiny':
        hicmat.matrix[hicmat.matrix == 0] = tiny
    if method == 'small':
        hicmat.matrix[hicmat.matrix == 0] = small
