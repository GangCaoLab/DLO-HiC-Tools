import matplotlib.pyplot as plt
import numpy as np

from dlo_hic.utils.wrap.hic import StrawWrap, CoolerWrap


MAIN_CHROM = ['chr'+str(i) for i in range(1, 23)] + ['chrX', 'chrY']


def matrix_val_range(arr, min_val=None, max_val=None):
    small = 1e-4
    arr_no_nan = arr[np.logical_not(np.isnan(arr))]

    if min_val is None:
        lt_min = arr[arr > arr.min()]
        if lt_min.shape[0] > 0:
            min_ = lt_min.min()
        else:
            min_ = small
    else:
        min_ = min_val

    if max_val is None:
        max_ = arr_no_nan.max()
    else:
        max_ = max_val

    if max_ <= min_:
        max_ = min_ + small

    return min_, max_



def plot_mat(mat, extent=None, color=None, log=True, val_range=(None, None)):
    fig, ax = plt.subplots()

    if color is None:  # default cmap
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('interaction',
                                                 ['#FFFFFF','#FFDFDF','#FF7575','#FF2626','#F70000'])
    else:
        cmap = plt.get_cmap(color)

    cmap.set_bad("white")
    cmap.set_under("white")

    c_min, c_max = matrix_val_range(mat, *val_range)

    img = ax.matshow(mat, cmap=cmap,
                     extent=extent,
                     aspect='auto')
    ax.xaxis.set_ticklabels([])

    import matplotlib.colors as colors
    if log:
        img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
    else:
        img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

    # plot colorbar
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size=0.12, pad=0.15)

    if log:
        from matplotlib.ticker import LogFormatter
        formatter = LogFormatter(10, labelOnlyBase=False)
        aa = np.array([1, 2, 5])

        def abs_inc(num):
            if num != 0:
                sign = num / abs(num)
                return int(sign * abs(num + 1))
            else:
                return 1

        lower_ = int(np.log10(c_min))
        upper_ = abs_inc(int(np.log10(c_max)))
        tick_values = np.concatenate([aa * 10 ** x for x in range(lower_, upper_)])

        c_bar = plt.colorbar(img, cax=cax, ticks=tick_values, format=formatter, fraction=0.98)
    else:
        c_bar = plt.colorbar(img, cax=cax, fraction=0.98)

    c_bar.solids.set_edgecolor("face")
    c_bar.ax.tick_params(labelsize='smaller')

    c_bar.ax.yaxis.set_ticks_position('right')
    return fig, ax



def plot_hic_mat(path, region, binsize='auto', balance=True):
    """
    plot a contact map of a chromosome.

    Args
    ----
    path : str
        Path to contact map file, in '.mcool' or '.hic' format.
    region : str
        Chromosome region, like: "chr1:1000-2000"
    binsize : int
        Resolution, use 'auto' to inference resolution automatically.
    balance : bool
        Use balanced matrix or not.
    color : str
        Color map.
    log : bool
        Use LogNorm or not.
    val_range : (float, float)
        Matrix value range.
    """
    from dlo_hic.utils.wrap.hic import GenomeRange
    wrap_func = StrawWrap if path.endswith(".hic") else CoolerWrap
    hic = wrap_func(path, binsize='auto', balance=balance)

    g_range = GenomeRange(region)
    mat = hic.fetch(g_range)
    fetched_binsize = hic.fetched_binsize
    #import ipdb; ipdb.set_trace()

    start, end = g_range.start, g_range.end
    extent = (start, end, end, start)
    fig, ax = plot_mat(mat, extent=extent, log=True)

    title = "region: \"" + region + "\"\nbinsize: " + str(fetched_binsize)
    ax.set_title(title)
    plt.tight_layout()

    return fig, ax


def plot_global_map(path, binsize='auto', balance=False, chroms=MAIN_CHROM):
    wrap_func = StrawWrap if path.endswith(".hic") else CoolerWrap
    hic = wrap_func(path, balance=balance)
    mat = hic.fetch_all(chroms=chroms, binsize=binsize)
#    mat[np.eye(mat.shape[0])==1] = 0

    max_ = mat.max() * 0.5

    fig, ax = plot_mat(mat, log=True, val_range=(None, max_))

    return fig, ax
    
