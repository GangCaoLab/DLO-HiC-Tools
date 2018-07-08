import os
from os.path import join
import sys
import logging
import multiprocessing

from .config import LOGGING_FMT, LOGGING_DATE_FMT

ncpu = multiprocessing.cpu_count()

class DLO_HiC(object):
    """
    Interface to DLO-HiC tools.

    Parameters
    ----------
    prefix : str
        The identifier of data and results. For example,
        your fastq data names "a1_R1.fastq",
        the prefix should be "a1". At the same time,
        the intermediate and final result's file preifx will also be "a1".

    workdir : str
        Path to working directory,
        final result and intermediate result will storaged in here.

    log_file : str, optional
        Path to log file, use '-' to stderr[default]

    log_level : int, optional
        Log level, default 10(DEBUG)

    processes : int, optional
        Use how many processes to run tools,
        defalut the cpu number in system.
    """

    DEFAULT_MISMATCH = 4
    DEFAULT_PHRED = 33
    DEFAULT_PET_LEN = 0

    DEFAULT_MAPQ = 1

    DEFAULT_THRESH_NUM_REST = 1
    DEFAULT_THRESH_SPAN = 50

    DEFAULT_THRESH_RR_DISTANCE = 0

    KEEP_RAW_PAIRS = False

    def __init__(self, prefix, workdir, log_file='-', log_level=10, processes=ncpu):
        self.prefix = prefix
        self.workdir = workdir
        self.processes = processes

        # set logger
        self.log = logging.getLogger("DLO_HiC<%s/%s>"%(workdir, prefix))
        if log_file == '-':
            handler = logging.StreamHandler(sys.stderr)
        else:
            handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
        self.log.addHandler(handler)
        self.log.setLevel(log_level)

        self.log.info("prefix(sample identifier): {}".format(self.prefix))
        self.log.info("working directory: {}".format(workdir))
        self.log.info("processes: {}".format(processes))

    def __check_required(self, para_name, kwargs):
        attr_name = para_name
        if not hasattr(self, attr_name):
            if para_name not in kwargs:
                msg = "parameter \"{}\" is need, you must provide it as an argument, or set it to attribute.".format(para_name)
                raise ValueError(msg)
            else:
                setattr(self, attr_name, kwargs[para_name])

    def __check_optional(self, para_name, kwargs, default):
        attr_name = para_name
        if not hasattr(self, attr_name):
            if para_name not in kwargs:
                msg = "parameter \"{}\" is not specified, set to default value: {}".format(para_name, default)
                self.log.info(msg)
                setattr(self, attr_name, default)
            else:
                setattr(self, attr_name, kwargs[para_name])

    def __mkdir(self, dirname):
        dir_ = join(self.workdir, dirname)
        if not os.path.exists(dir_):
            os.mkdir(dir_)
        return dir_

    def extract_PET(self, *args, **kwargs):
        from dlo_hic.tools.main_process.extract_PET import main

        self.__check_required('fastq', kwargs)
        self.__check_required('linker_a', kwargs)
        self.__check_optional('linker_b', kwargs, None)
        self.__check_optional('mismatch', kwargs, DLO_HiC.DEFAULT_MISMATCH)
        self.__check_required('rest', kwargs)
        self.__check_optional('phred', kwargs, DLO_HiC.DEFAULT_PHRED)
        self.__check_optional('PET_len', kwargs, DLO_HiC.DEFAULT_PET_LEN)

        dir_ = self.__mkdir('01-extpet')

        PET1_file = join(dir_, self.prefix+'.pet1.fq')
        PET2_file = join(dir_, self.prefix+'.pet2.fq')
        PET_log = join(dir_, self.prefix+'.pet.log')

        main(self.fastq, PET1_file, PET2_file,
             self.linker_a, self.linker_b,
             self.mismatch, self.rest, self.phred, self.processes, self.PET_len,
             PET_log)

        setattr(self, 'PET1_file', PET1_file)
        setattr(self, 'PET2_file', PET2_file)
        setattr(self, 'PET_log', PET_log)

    def build_bedpe(self, *args, **kwargs):
        from dlo_hic.tools.main_process.build_bedpe import main

        self.__check_required('PET1_file', kwargs)
        self.__check_required('PET2_file', kwargs)
        self.__check_required('bwa_index', kwargs)
        self.__check_optional('mapq', kwargs, DLO_HiC.DEFAULT_MAPQ)

        dir_ = self.__mkdir('02-bedpe')

        bedpe_file = join(dir_, self.prefix+'.uniq.bedpe')
        bwa_log_file = join(dir_, self.prefix+'.bwa.log')

        main('fastq', self.PET1_file, self.PET2_file, bedpe_file,
             self.processes, self.bwa_index, self.mapq, bwa_log_file)
        
        setattr(self, 'bedpe_file', bedpe_file)
        setattr(self, 'bwa_log_file', bwa_log_file)

    def noise_reduce(self, *args, **kwargs):
        from dlo_hic.tools.main_process.noise_reduce import main

        self.__check_required('bedpe_file', kwargs)
        self.__check_required('restriction_site_file', kwargs)
        self.__check_optional('threshold_num_rest', kwargs, DLO_HiC.DEFAULT_THRESH_NUM_REST)
        self.__check_optional('threshold_span', kwargs, DLO_HiC.DEFAULT_THRESH_SPAN)

        dir_ = self.__mkdir('03-nr')
        nr_file = join(dir_, self.prefix+'.nr.bedpe')

        main(self.bedpe_file, nr_file, 
             self.restriction_site_file, self.processes,
             self.threshold_num_rest, self.threshold_span)

        setattr(self, 'nr_file', nr_file)

    def remove_reduncancy(self, *args, **kwargs):
        from dlo_hic.tools.main_process.remove_redundancy import main

        self.__check_required('nr_file', kwargs)
        self.__check_optional('threshold_rr_distance', kwargs, DLO_HiC.DEFAULT_THRESH_RR_DISTANCE)
        dir_ = self.__mkdir('04-rr')
        rr_file = join(dir_, self.prefix+'.rr.bedpe')
        main(self.nr_file, rr_file, self.threshold_rr_distance)
        setattr(self, 'rr_file', rr_file)

    def bedpe2pairs(self, *args, **kwargs):
        from dlo_hic.tools.main_process.bedpe2pairs import main

        self.__check_required('rr_file', kwargs)
        self.__check_optional('keep_raw_pairs', kwargs, DLO_HiC.KEEP_RAW_PAIRS)

        dir_ = self.__mkdir('05-rr')
        pairs_file = join(dir_, self.prefix+".pairs")

        main(self.rr_file, pairs_file, self.keep_raw_pairs)

        setattr(self, 'pairs_file', pairs_file)
