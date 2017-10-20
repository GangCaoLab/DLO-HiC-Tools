import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage


def plot_TAD(chrmat, chr_, start=None, end=None,
        transform=np.log10, cmap='RdYlBu_r', cbar=True, figsize=(10, 10), origin='lower'):
    """
    plot TAD.
    :chrmat: HicChrMatrix object
    :chr_: TAD's chromosome
    :start: TAD's start position
    :end: TAD's end position

    if start or end is 'False' or 'None', will plot whole chromosome.

    """
    mat = chrmat.chromosome(chr_).matrix

    if transform:
        with np.errstate(divide='ignore'):
            mat = transform(mat)
        mat[mat==-np.inf] = 0

    if start and end:
        # cut matrix to range: start <-> end
        if not start < end:
            raise ValueError("TAD start position must less than end position.")
        start_, end_ = start//chrmat.bin_size, end//chrmat.bin_size # calculate offset
        mat = mat[start_:end_+1, start_:end_+1]
    else:
        start = 0
        length = dict(chrmat.chr_len)[chr_]
        end = length - 1

    def rotate_and_cut(mat):
        res = ndimage.rotate(mat, 45)
        h = res.shape[0]
        res = res[h//2:, :]
        return res

    tad = rotate_and_cut(mat)

    # remove points outside TAD
    mask = np.ones(mat.shape)
    mask = rotate_and_cut(mask)
    tad = np.where(mask == 0, np.nan, tad)

    fig, ax = plt.subplots(figsize=figsize)

    img = ax.imshow(tad, cmap=cmap, origin=origin)

    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(img, cax=cax)

    ax.set_xticks([0, tad.shape[1]])
    ax.set_yticks([])
    ax.set_xticklabels([str(start), str(end)])
    
    # hide spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    return img
