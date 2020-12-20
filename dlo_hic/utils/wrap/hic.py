import logging
log = logging.getLogger(__name__)

from abc import ABC, abstractmethod


MAIN_CHROM = ['chr'+str(i) for i in range(1, 23)] + ['chrX']


class HiCWrap(ABC):


    @abstractmethod
    def fetch(self, genome_range1, genome_range2, binsize):
        pass

    @property
    @abstractmethod
    def chromosomes(self):
        pass

    @property
    @abstractmethod
    def chrom_length(self):
        pass

    @property
    @abstractmethod
    def resolutions(self):
        pass

    def fetch_chrom(self, chrom, binsize='auto'):
        length = self.chrom_length[chrom]
        g_range = chrom + ":1" + "-" + str(length)
        mat = self.fetch(g_range, binsize=binsize)
        return mat, length

    def fetch_all(self, chroms=MAIN_CHROM, binsize='auto', bin_num=2000):
        total_length = sum([v for k,v in self.chrom_length.items()])
        if binsize == 'auto':
            binsize = infer_resolution(total_length, self.resolutions, bin_num)

        chrom_length = [(k, v) for k,v in self.chrom_length.items() if k in chroms]

        import numpy as np
        cache = {}
        rows = []
        for chrom_1, length_1 in chrom_length:
            row = []
            for chrom_2, length_2 in chrom_length:
                g_range_1 = chrom_1 + ":1" + "-" + str(length_1)
                g_range_2 = chrom_2 + ":1" + "-" + str(length_2)
                if (g_range_2, g_range_1) in cache:
                    mat = cache[(g_range_2, g_range_1)].T
                else:
                    mat = self.fetch(g_range_1, g_range_2, binsize=binsize)
                    cache[(g_range_1, g_range_2)] = mat
                row.append(mat)
            rows.append(np.hstack(row))
        res = np.vstack(rows)

        mat_widths = [c.shape[1] for c in row]

        chrom_start_pos = []
        start = 0
        for (chrom, length), width in zip(chrom_length, mat_widths):
            chrom_start_pos.append((chrom, start))
            start += width

        return res, chrom_start_pos, binsize


class StrawWrap(HiCWrap):
    """
    A wrap for straw Python API, for read .hic file.

    Parameters
    ----------
    path : str
        Path to the '.hic' file

    balance : {bool, 'VC', 'VC_SQRT', 'KR'}
        Method for the matrix normalization.
        default 'KR'

    binsize : {'auto', int}
        resolution of the data. for example 5000.
        'auto' for calculate resolution automatically.
        default 'auto'

    """
    def __init__(self, path, balance='KR', binsize='auto'):
        self.hic_file = path
        if balance is True:
            normalization = 'KR'  # default normalization method
        elif balance is False:
            normalization = 'NONE'
        else:
            normalization = balance
        self.normalization = normalization
        self.binsize = binsize
        self._chromosomes, self._resolutions, self.masterindex, self.genome, self.metadata = self.__info()
        self.fetched_binsize = None

    @property
    def chromosomes(self):
        return [v[1] for _,v in self._chromosomes.items()]

    @property
    def resolutions(self):
        return sorted(self._resolutions, reverse=True)

    @property
    def chrom_length(self):
        from collections import OrderedDict
        res = OrderedDict()
        for k, v in self._chromosomes.items():
            if k != 0:
                chr_ = v[1] if v[1].startswith('chr') else change_chrom_names(v[1])
                res[chr_] = v[2]
        return res

    def fetch(self, genome_range1, genome_range2=None, binsize=None):
        """
        Return
        ------
        matrix : numpy.ndarray
        """

        if genome_range2 is None:
            genome_range2 = genome_range1

        if isinstance(genome_range1, str):
            genome_range1 = GenomeRange(genome_range1)
        if isinstance(genome_range2, str):
            genome_range2 = GenomeRange(genome_range2)

        if genome_range1.chrom not in self.chromosomes:
            genome_range1.change_chrom_names()
        if genome_range2.chrom not in self.chromosomes:
            genome_range2.change_chrom_names()

        if binsize is None:
            binsize = self.__infer_binsize(genome_range1)
            self.fetched_binsize = binsize  # expose fetched binsize

        straw_list = self.__fetch_straw_list(genome_range1, genome_range2, binsize)
        matrix = self.__list_to_matrix(straw_list, genome_range1, genome_range2, binsize)
        return matrix

    def __infer_binsize(self, genome_range):
        if self.binsize == 'auto':
            binsize = infer_resolution(genome_range.length, self.resolutions)
        else:
            binsize = self.binsize
        return binsize

    def __fetch_straw_list(self, genome_range1, genome_range2, binsize):
        from .straw import straw
        chr1loc = str(genome_range1).replace('-', ':')
        chr2loc = str(genome_range2).replace('-', ':')
        slist = straw(self.normalization, self.hic_file, chr1loc, chr2loc, 'BP', binsize)
        return slist

    def __list_to_matrix(self, straw_list, genome_range1, genome_range2, binsize):
        import numpy as np

        #flag = False
        #if genome_range2.length > genome_range1.length:
        #    flag = True
        #    genome_range1, genome_range2 = genome_range2, genome_range1

        binlen1 = (genome_range1.length // binsize) + 1
        binlen2 = (genome_range2.length // binsize) + 1

        mat = np.zeros((binlen1, binlen2), dtype=np.float64)
        for loc1, loc2, c in zip(*straw_list):
            bin1id = min((loc1 - genome_range1.start) // binsize, binlen1 - 1)
            bin2id = min((loc2 - genome_range2.start) // binsize, binlen2 - 1)
            if binlen1 == binlen2:
                mat[bin1id, bin2id] = c
                mat[bin2id, bin1id] = c
            else:
                mat[bin1id, bin2id] = c

        #if flag:
        #    mat = mat.T

        return mat

    def __info(self):
        """
        from hic2cool code:
            https://github.com/4dn-dcic/hic2cool/blob/master/hic2cool/hic2cool_utils.py#L73

        Takes in a .hic file and returns a dictionary containing information about
        the chromosome. Keys are chromosome index numbers (0 through # of chroms
        contained in file) and values are [chr idx (int), chr name (str), chrom
        length (str)]. Returns the masterindex used by the file as well as the open
        file object.
        """
        import sys
        import struct
        from .straw import readcstr

        req = open(self.hic_file, 'rb')

        chrs = {}
        resolutions = []
        magic_string = struct.unpack(b'<3s', req.read(3))[0]
        req.read(1)
        if (magic_string != b"HIC"):
            print('This does not appear to be a HiC file; '
                  'magic string is incorrect')
            sys.exit()
        global version
        version = struct.unpack(b'<i', req.read(4))[0]
        masterindex = struct.unpack(b'<q', req.read(8))[0]
        genome = b""
        c = req.read(1)
        while (c != b'\0'):
            genome += c
            c = req.read(1)
        genome = genome.decode('ascii')
        # metadata extraction
        metadata = {}
        nattributes = struct.unpack(b'<i', req.read(4))[0]
        for x in range(nattributes):
            key = readcstr(req)
            value = readcstr(req)
            metadata[key] = value
        nChrs = struct.unpack(b'<i', req.read(4))[0]
        for i in range(0, nChrs):
            name = readcstr(req)
            length = struct.unpack(b'<i', req.read(4))[0]
            if name and length:
                chrs[i] = [i, name, length]
        nBpRes = struct.unpack(b'<i', req.read(4))[0]
        # find bp delimited resolutions supported by the hic file
        for x in range(0, nBpRes):
            res = struct.unpack(b'<i', req.read(4))[0]
            resolutions.append(res)
        return chrs, resolutions, masterindex, genome, metadata


class CoolerWrap(HiCWrap):
    """
    wrap for cooler file,
    deal with multi resolution.

    Parameters
    ----------
    path : str
        Path to the cooler file.

    binsize : {'auto', int}
        resolution of the data. for example 5000.
        'auto' for calculate resolution automatically.
        default 'auto'

    balance : bool
        Balance the matrix or not.
        default True
    """
    def __init__(self, path, binsize='auto', balance=True):
        import cooler
        self.path = path

        self.is_multi = is_multi_cool(path)
        if self.is_multi:
            self.coolers = self.__load_multi_coolers(path)
        else:
            self.cool = cooler.Cooler(path)

        self.binsize = binsize
        self.balance = balance

    @property
    def chromosomes(self):
        if self.is_multi:
            c = self.coolers[self.resolutions[0]]
        else:
            c = self.cool
        return [k for k,_ in dict(c.chromsizes).items()]

    @property
    def chrom_length(self):
        if self.is_multi:
            c = self.coolers[self.resolutions[0]]
        else:
            c = self.cool
        from collections import OrderedDict
        res =  OrderedDict()
        for k, v in c.chromsizes.iteritems():
            chr_ = k if k.startswith('chr') else change_chrom_names(k)
            res[chr_] = v
        return res

    @property
    def resolutions(self):
        resolutions = get_cooler_resolutions(self.path, self.is_multi)
        return sorted(resolutions, reverse=True)

    def __load_multi_coolers(self, path):
        from collections import OrderedDict
        import cooler

        def cooler_reso(resolution):
            from h5py import File
            with File(path, 'r') as f:
                if "resolutions" in f:
                    c = cooler.Cooler(path + "::/resolutions/{}".format(resolution))
                else:
                    for grp_name in f:
                        grp = f[grp_name]
                        if str(grp.attrs['bin-size']) == str(resolution):
                            c = cooler.Cooler(path + "::{}".format(grp_name))
                            break
            return c

        coolers = OrderedDict([(reso, cooler_reso(reso)) for reso in self.resolutions])
        return coolers

    def get_cool(self, genome_range):
        if self.is_multi:
            resolutions = [k for k in self.coolers.keys()]
            if self.binsize == 'auto':
                binsize = infer_resolution(genome_range.length, resolutions)
            else:
                assert self.binsize in resolutions, \
                    "Multi-Cooler file not contain the resolution {}.".format(self.binsize)
                binsize = int(self.binsize)
            self.fetched_binsize = binsize
            cool = self.coolers[binsize]
        else:
            cool = self.cool
            self.fetched_binsize = cool.binsize  # expose fetched binsize
        return cool

    def fetch(self, genome_range1, genome_range2=None, binsize=None):
        if genome_range2 is None:
            genome_range2 = genome_range1

        if isinstance(genome_range1, str):
            genome_range1 = GenomeRange(genome_range1)
        if isinstance(genome_range2, str):
            genome_range2 = GenomeRange(genome_range2)

        if genome_range1.chrom not in self.chromosomes:
            genome_range1.change_chrom_names()
        if genome_range2.chrom not in self.chromosomes:
            genome_range2.change_chrom_names()

        if binsize is not None:
            if self.is_multi:
                cool = self.coolers[binsize]
            else:
                cool = self.cool
        else:
            cool = self.get_cool(genome_range1)

        try:
            mat = cool.matrix(balance=self.balance).fetch(str(genome_range1), str(genome_range2))
        except ValueError as e:
            log.warning(str(e))
            log.warning("Data is not balanced, force to use unbalanced matrix.")
            mat = cool.matrix(balance=False).fetch(str(genome_range1), str(genome_range2))

        return mat


def file_type(path):
    if path.endswith(".hic"):
        return '.hic'
    else:
        p = path.split("::")[0]
        if p.endswith(".cool") or p.endswith(".mcool"):
            return '.cool'
        else:
            raise NotImplementedError("file type of {} not supported yet".format(path))


def infer_resolution(range_length, resolutions, bin_thresh=500):
    """
    Inference appropriate resolution.
    """
    resolutions = sorted(resolutions)
    reso = resolutions[0]
    for r in resolutions:
        num_bins = range_length // r
        if num_bins >= bin_thresh:
            reso = r
        else:
            break
    return reso


def is_multi_cool(cooler_file):
    """
    Judge a cooler is muliti-resolution cool or not.
    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import re
    if re.match(".+::.+$", cooler_file):
        return False

    import h5py
    h5_file = h5py.File(cooler_file, 'r')
    is_multi = 'pixels' not in h5_file  # use "pixels" group distinguish is multi-cool or not
    h5_file.close()
    return is_multi


def get_cooler_resolutions(cooler_file, is_multi=True):
    """
    Get the resolutions of a muliti-cooler file
    Parameters
    ----------
    cooler_file : str
        Path to cooler file.
    """
    import h5py
    h5_file = h5py.File(cooler_file, 'r')
    if is_multi:
        if 'resolutions' in h5_file:
            resolutions = list(h5_file['resolutions'])
            resolutions = [int(res) for res in resolutions]
        else:
            resolutions = [int(h5_file[i].attrs['bin-size']) for i in list(h5_file)]
        resolutions.sort()
        h5_file.close()
    else:
        resolutions = [int(h5_file.attrs['bin-size'])]

    return resolutions


class GenomeRange(object):
    """
    Express a range on the genome.
    Attributes
    ----------
    chrom : str
        chromosome
    start : int
        start position
    end : int
        end position
    """

    def __init__(self, *args):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> str(range1)
        'chr1:1000-2000'
        >>> range2 = GenomeRange("chr2:2000-4000")
        >>> (range2.chrom, range2.start, range2.end)
        ('chr2', 2000, 4000)
        >>> range3 = GenomeRange("chr1", 2000, 1000)
        Traceback (most recent call last):
        ...
        ValueError: Please check that the region end is larger than the region start. Values given: start: 2000, end: 1000
        """
        if len(args) == 1:
            chrom, start, end = GenomeRange.parse_region_string(args[0])
        elif len(args) == 3:
            chrom, start, end = args
        else:
            raise ValueError("inappropriate init arguments. "
                             "correct example: `range1 = GenomeRange(\"chr1:1000-2000\")` or "
                             "`range1 = GenomeRange(\"chr1\", 1000, 2000)`")

        if end < start:
            raise ValueError("Please check that the region end is larger than the region start. "
                             "Values given: start: {}, end: {}".format(start, end))

        self.chrom = chrom
        self.start = start
        self.end = end

    @staticmethod
    def parse_region_string(region_string):
        """
        splits a region string into
        a (chrom, start, end) tuple
        Parameters
        ----------
        region_string : str
            Region string to be parsed, like: "chr:start-end"
        Return
        ------
        result : tuple of str
            Result tuple (chrom, start, end)
        >>> GenomeRange.parse_region_string("chr1:10-20")
        ('chr1', 10, 20)
        >>> GenomeRange.parse_region_string("chr1:0")
        Traceback (innermost last):
         ...
        ValueError: Failure to parse region string, please check that region format should be like "chr:start-end".
        """
        if region_string:
            # separate the chromosome name and the location using the ':' character
            chrom, position = region_string.strip().split(":")

            # clean up the position
            for char in ",.;|!{}()":
                position = position.replace(char, '')

            position_list = position.split("-")
            try:
                region_start = int(position_list[0])
                region_end = int(position_list[1])
            except:
                raise ValueError("Failure to parse region string, please check that region format "
                                 "should be like \"chr:start-end\".")

            return chrom, region_start, region_end

    def change_chrom_names(self):
        """
        >>> range1 = GenomeRange("chr1", 1000, 2000)
        >>> range1.chrom
        'chr1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        '1'
        >>> range1.change_chrom_names()
        >>> range1.chrom
        'chr1'
        """
        self.chrom = change_chrom_names(self.chrom)

    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)

    @property
    def length(self):
        """
        >>> range1 = GenomeRange("chr1", 0, 1000)
        >>> range1.length
        1000
        """
        return self.end - self.start

    def __eq__(self, other):
        """
        >>> GenomeRange('chr1', 1000, 2000) == GenomeRange("chr1:1000-2000")
        True
        >>> GenomeRange("chr1:1000-2000") == GenomeRange("1:1000-2000")
        False
        """
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def __contains__(self, another):
        if another.chrom != self.chrom:
            return False
        if another.start < self.start:
            return False
        if another.end > self.end:
            return False
        return True


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    >>> change_chrom_names("chr1")
    '1'
    >>> change_chrom_names("1")
    'chr1'
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom
