# -*- coding: utf-8 -*-


import os
import sys
import cPickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from .plot.matrix import plot_hicmat, plot_chrmat


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

    def save_matrix(self, filename):
        """ save numpy matrix/ndarray to file. """
        np.save(filename, self.matrix)

    @staticmethod
    def save(hicmat, fname):
        """ 
        save object to two files:
            .npy for store numpy matrix
            .pyobj file for store class instance
        using numpy and cPickle. 

        >>> HicMatrix.save(hicmat, "test")

        $ ls
        test.npy test.pyobj

        """
        prefix = os.path.splitext(fname)[0]
        hicmat.save_matrix(prefix)
        with open(prefix + ".pyobj", 'w') as f:
            mat = hicmat.matrix
            hicmat.matrix = None
            cPickle.dump(hicmat, f)
            hicmat.matrix = mat

    @staticmethod
    def load(fname):
        """ 
        load object from .npy file and .pyobj file.

        >>> HicMatrix.load("test")
        array([[1 2 0 ... 1]
                ...
               [1 3 0 ... 1]])
        
        """
        prefix = os.path.splitext(fname)[0]
        try:
            with open(prefix + ".pyobj") as f:
                hicmat = cPickle.load(f)
            hicmat.matrix = np.load(prefix + ".npy")
        except IOError as ioe:
            sys.stderr.write("[warning] Unable to open .pyobj file " + str(ioe) + " ")
            sys.stderr.write("will creat HicMatrix recover only from numpy matrix file.\n")
            hicmat = HicMatrix(np.load(prefix + ".npy"))
        return hicmat

    def to_file(self, filename):
        """ 
        save object to binary file. 
        equal to HicMatrix.save

        >>> hicmat.save("test")
        
        """
        HicMatrix.save(self, filename)

    def __repr__(self):
        return "HicMatrix: \n" + repr(self.matrix)


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
            can use HicChrMatrix.load_chr_len load from tab split file.

        """
        self.bin_size = bin_size

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

    @staticmethod
    def load_chr_len(file_chr_len):
        """ 
        load a list of (chromesome, length) pair.

        :file_chr_len: a Tab split file record,
        the order and length(in basepair) of chrosome, file format like this:

        chr1    248956422
        chr2    242193529
        ...
        chrM    16569
        
        """
        chr_len = list()
        with file_chr_len as f:
            for line in f:
               chr_, len_ = line.strip().split("\t") 
               len_ = int(len_)
               chr_len.append((chr_, len_))
        return chr_len

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

    def locate(self, chr_x, x, chr_y, y):
        """
        locate an interaction on matrix.

        firstly, find bin_x and bin_y, then
        interaction strength between bin_x, bin_y increase one.

        matrix[bin_x, bin_y] += 1
        matrix[bin_y, bin_x] += 1
        
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
        self.matrix[abs_x, abs_y] += 1
        self.matrix[abs_y, abs_x] += 1

    def plot(self, *args, **kwargs):
        """ plot matrix with chromosomes information. """
        img = plot_chrmat(self, *args, **kwargs)
        return img
