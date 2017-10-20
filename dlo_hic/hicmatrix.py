# -*- coding: utf-8 -*-

import copy
from functools import wraps

import numpy as np

from .plot.matrix import plot_hicmat, plot_chrmat
from .IO.matrix import save_hicmat


def mat_operation(func):
    """ Decorator for HicMatrix operation method. """
    @wraps(func)
    def warp(self, other):
        if self.matrix.shape != other.matrix.shape:
            raise ValueError("Two matrix must in same shape.")
        res = copy.copy(self)
        res.matrix = func(self.matrix, other.matrix)
        return res
    return warp


class HicMatrix:
    """
    The basic abstrction of Hic result matrix.
    provide some basic operation (e.g. plot, write and load)
    for hic result matrix.
    
    """

    def __init__(self, matrix):
        """
        :matrix: numpy matrix or ndarray, dtype is numpy.int

        """
        self.matrix = matrix
 
    def plot(self, *args, **kwargs):
        """ 
        plot the matrix.
        :transform: scale transformation function, default False.
            if you want to transform scale, can specify some numpy function e.g. np.log2 .

        """
        img = plot_hicmat(self, *args, **kwargs)
        return img

    def save(self, *args, **kwargs):
        """ 
        save object to npz files:
        """
        save_hicmat(self, *args, **kwargs)

    def __repr__(self):
        return "HicMatrix: \n" + repr(self.matrix)

    @mat_operation
    def __add__(a, b):
        return a + b

    @mat_operation
    def __sub__(a, b):
        return a - b

    @mat_operation
    def __div__(a, b):
        return a / b

    @mat_operation
    def __mul__(a, b):
        return a * b

    @mat_operation
    def log2_fold_change(a, b):
        """
        get the log2(Fold Change) from self to other
        a.log2_fold_change(b) = log2(b / a)
        """
        fc = b / a
        log2_fc = np.log2(fc)
        return log2_fc


class HicChrMatrix(HicMatrix):
    """
    The extend of HicMatrix,
    contain the chromosomes information.
    
    This is more convenient for build interaction matrix.
    
    >>> hicmat = HicChrMatrix(chr_len, 10000)
    >>> hicmat
    array([0, 0, ..., 0, 0],
          ...,
          [0, 0, ..., 0, 0])
    >>> hicmat.locate(("chr1", 200, "chr1", 100))
    array([1, 0, ..., 0, 0],
          ...,
          [0, 0, ..., 0, 0])

    """
    def __init__(self, chr_len, bin_size):
        """
        :chr_len: a list of (chromosome_name, length) pair,
            record both order and length(in basepair) of chromosomes.

        :bin_size: the length of bins in the matrix.
            can use dlo_hic.IO.load_chr_len load from tab split file.

        """
        self.bin_size = bin_size
        self.chr_len = chr_len

        # build a linner space to represent all chromosomes bin range
        chromosomes, lengths = [], []
        axis, pos = dict(), 0
        for chr_, len_ in chr_len:
            num_bins = len_ // bin_size
            if num_bins == 0: # chromosome length is too short
                continue
            else:
                axis[chr_] = (pos, pos+num_bins-1)
                pos += num_bins
                chromosomes.append(chr_)
                lengths.append(num_bins)
        self.axis = axis
        self.chromosomes = chromosomes
        self.lengths = np.asarray(lengths)

        # number of all bins
        num_bins = sum([len_//bin_size for _, len_ in chr_len])
        self.num_bins = num_bins
        # init super class
        self.matrix = np.zeros([num_bins, num_bins], dtype=np.int)
        HicMatrix.__init__(self, self.matrix)

    def chromosome(self, chr_):
        """
        return a sub chromosome's hic matrix.

        >>> chr1_mat = hicmat.chromosome("chr1")

        """
        chr_span = self.axis[chr_]
        chr_matrix = self.matrix[chr_span[0]:chr_span[1]+1,
                                 chr_span[0]:chr_span[1]+1]
        return HicMatrix(chr_matrix)

    def __getitem__(self, (chr_a, chr_b)):
        """
        a subparts of HicChrMatrix: chr_a and chr_b interaction.

        >>> chr1_chr2_mat = hicmat["chr1", "chr2"]

        """ 
        chr_a_span = self.axis[chr_a]
        chr_b_span = self.axis[chr_b]
        ab_matrix = self.matrix[chr_a_span[0]:chr_a_span[1]+1,
                                chr_b_span[0]:chr_b_span[1]+1]
        return HicMatrix(ab_matrix)

    def __repr__(self):
        return "HicChrMatrix: \n" + repr(self.matrix)

    def locate(self, chr_x, x, chr_y, y, val=1):
        """
        locate an interaction on matrix.

        firstly, find bin_x and bin_y, then
        interaction strength between bin_x, bin_y increase val(default 1).

        matrix[bin_x, bin_y] += val
        matrix[bin_y, bin_x] += val
        
        >>> hicmat.bin_size
        10
        >>> hicmat.matrix
        array([[1, 2]
               [2, 0]])
        >>> hicmat.locate(5, 17)
        >>> hicmat.matrix
        array([[1, 3]
               [3, 0]])
        """
        chr_span_x = self.axis[chr_x]
        chr_span_y = self.axis[chr_y]
        # the offset of x and y in matrix
        offset_x = x//self.bin_size
        offset_y = y//self.bin_size
        # the absulute adderss of x and y in matrix
        abs_x = offset_x + chr_span_x[0]
        abs_y = offset_y + chr_span_y[0]
        self.matrix[abs_x, abs_y] += val
        self.matrix[abs_y, abs_x] += val

    def plot(self, *args, **kwargs):
        """ plot matrix with chromosomes information. """
        img = plot_chrmat(self, *args, **kwargs)
        return img
