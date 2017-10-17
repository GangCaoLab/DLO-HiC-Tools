import matplotlib.pyplot as plt 
from matplotlib import colors

def plot_hicmat(hicmat, transform=False, cmap="RdBu_r", cbar=False, figsize=(10, 10),
        norm=colors.SymLogNorm(1)):
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
    img = ax.imshow(mat, origin="lower", cmap=cmap,
            extent=(0, mat.shape[0], 0, mat.shape[1]),
            norm=norm)
    if cbar:
        cb = fig.colorbar(img)
    return img


def plot_chrmat(chrmat, transform=False, cmap="RdBu_r", cbar=True,
        figsize=(20, 14), ticks=True, norm=colors.SymLogNorm(1)):
    """ plot matrix with chromosomes information. """
    img = HicMatrix.plot(self, transform=transform, cmap=cmap, cbar=cbar,\
            figsize=figsize, norm=norm)
    ax = img.axes
    for i in self.lengths.cumsum(): # plot line between chromosomes
        ax.axhline(i, linewidth=1, color="#222222")
        ax.axvline(i, linewidth=1, color="#222222")
    if ticks:
        ticks = np.append(np.zeros([1]), self.lengths.cumsum())
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(self.chromosomes, rotation="vertical")
        ax.set_yticklabels(self.chromosomes)
    return img
