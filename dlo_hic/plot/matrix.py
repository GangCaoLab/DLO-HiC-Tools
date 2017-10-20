import numpy as np
import matplotlib.pyplot as plt 


def plot_hicmat(hicmat, transform=np.log10, cmap="RdYlBu_r",
        cbar=False, figsize=(10, 10)):
    """ 
    plot the matrix.
    :transform: scale transformation function, default False.
        if you want to transform scale, can specify some numpy function e.g. np.log2 .

    """
    if transform:
        with np.errstate(divide='ignore'):
            mat = transform(hicmat.matrix)
        mat[mat==-np.inf] = 0
    else:
        mat = hicmat.matrix
    fig, ax = plt.subplots(figsize=figsize)
    img = ax.imshow(mat, origin="lower",
            cmap=cmap, extent=(0, mat.shape[0], 0, mat.shape[1]))
    #img = ax.imshow(mat, origin="lower", cmap=cmap,
    #        extent=(0, mat.shape[0], 0, mat.shape[1]),
    #        norm=norm)
    if cbar:
        cb = fig.colorbar(img)
    return img


def plot_chrmat(chrmat, transform=np.log10, cmap="RdYlBu_r", cbar=True,
        figsize=(20, 14), ticks=True):
    """ plot matrix with chromosomes information. """
    img = plot_hicmat(chrmat, transform=transform, cmap=cmap, cbar=cbar,
            figsize=figsize)
    ax = img.axes
    for i in chrmat.lengths.cumsum(): # plot line between chromosomes
        ax.axhline(i, linewidth=1, color="#222222", alpha=0.5)
        ax.axvline(i, linewidth=1, color="#222222", alpha=0.5)
    if ticks:
        xticks = np.append(np.zeros([1]), chrmat.lengths.cumsum())
        yticks = xticks[::-1]
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels(chrmat.chromosomes, rotation="vertical")
        ax.set_yticklabels([""] + chrmat.chromosomes[::-1])
    return img
