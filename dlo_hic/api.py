import sys
import logging
import multiprocessing

from .config import LOGGING_FMT, LOGGING_DATE_FMT

ncpu = multiprocessing.cpu_count()

class DLO_HiC(object):
    """
    Interface to DLO-HiC tools.
    """
    def __init__(self, prefix, workdir, log_file='-', log_level=10, processes=ncpu):
        """
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
            Path to log file, use '-' to stdout[default]
        log_level : int, optional
            Log level, default 10(DEBUG)
        processes : int, optional
            Use how many processes to run tools,
            defalut the cpu number in system.
        """
        self.prefix = prefix
        self.workdir = workdir
        self.processes

        # set logger
        self.log = logging.getLogger("DLO_HiC<%s/%s>"%(workdir, prefix))
        if log_file == '-':
            handler = logging.StreamHandler(sys.stderr)
        else:
            handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter(fmt=LOGGING_DATE_FMT, datefmt=LOGGING_DATE_FMT))
        self.log.addHandler(handler)
        self.log.setLevel(log_level)


    def extract_PET(self, fastq, linker_a, linker_b="None"):
        pass
        #from dlo_hic.tools.main_process.extract_PET import main

        #main(fastq, PET1, PET2,
        #    linker_a, linker_b,
        #    mismatch, rest, phred, processes, PET_len)

        #self.PET1 = PET1
        #self.PET2 = PET2
        