import os
import sys
from datetime import datetime
import logging
import subprocess
from subprocess import Popen, PIPE
import multiprocessing

log = logging.getLogger(__name__)

cores = multiprocessing.cpu_count()

class BWA():
    """
    bwa command wrap.
    now, only use for single end alignment.

    Parameters
    ----------
    index_prefix : str
        path prefix to the bwa index.

    algorithm : {'mem', 'aln'}
        bwa algorithm.
    
    log_file : {str, '-'}
        path to log file, use '-' for output log to stdout.
    """
    def __init__(self, index_prefix, algorithm='aln', log_file='-'):
        self.index_prefix = index_prefix
        self.algorithm = algorithm
        if log_file == '-':
            self.log_file = sys.stdout
        else:
            self.log_file = open(log_file, 'a')

    def _check_call(self, cmd):
        """ 
        call cmd and log to self.log_file
        """
        process = Popen(cmd, stderr=PIPE, shell=True)
        with process.stderr:
            for line in process.stderr:
                line = line.decode("utf-8")
                timestr = datetime.now().strftime("%m/%d/%y %H:%M:%S")
                line = "[{}]".format(timestr) + line # add datetime information
                self.log_file.write(line)
        exitcode = process.wait() 
        if exitcode != 0:
            raise subprocess.CalledProcessError(exitcode, cmd)

    def _stream(self, cmd):
        """
        call cmd, yiels stdout, log to self.log_file
        """
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        with process.stdout as out, process.stderr as err:
            while True:
                log_line = err.readline()
                if log_line:
                    log_line = log_line.decode("utf-8")
                    timestr = datetime.now().strftime("%m/%d/%y %H:%M:%S")
                    log_line = "[{}]".format(timestr) + log_line # add datetime information
                    self.log_file.write(log_line)
                line = out.readline()
                if not line:
                    break
                line = line.decode("utf-8")
                yield line

    def index(self, ref_fasta):
        """ build index """
        index_prefix = self.index_prefix
        msg = "build bwa index. fasta file: %s, index_prefix: %s"%(ref_fasta, index_prefix)
        log.info(msg)
        cmd = "bwa index -p {} {}".format(index_prefix, ref_fasta)
        self._check_call(cmd)

    def run(self, input, output_prefix,
            thread=cores, bam=True,
            max_diff=None,
            params=None, samse_params=None):
        """
        launch bwa aln/mem command, return process

        Parameters
        ----------
        input : str
            Path to input fastq file.
        output_prefix : str
            Path prefix to output file.
        thread : int
            Number of cores used for run bwa.
        bam : bool
            convert output to bam or not.
        max_diff : {int, None}
            set "-n" parameter.
        params : {str, None}
            Additional bwa alignment params.
        samse_params : {str, None}
            Additional samse params.
        """
        # check index exist or not
        index_prefix = self.index_prefix
        suffix = ["amb", "bwt", "sa", "ann", "pac"]
        files = [index_prefix + "." + s for s in suffix]
        exists = [os.path.exists(f) for f in files]
        all_exists = all(exists)
        if not all_exists:
            raise IOError('bwa index not found, run BWA.index firstly.')

        mem = True if self.algorithm == "mem" else False
        if mem:  # bwa mem
            log.info("run bwa mem algorithm on %s"%input)
            cmd = "bwa mem {} {} -t {}".format(self.index_prefix, input, thread)
            if max_diff is not None:
                cmd += " -n {} ".format(max_diff)
            if params:
                cmd += " " + params
            if bam:
                cmd += " | samtools view -bh"

            # run cmd
            suffix = "sam" if not bam else "bam"
            self._check_call(cmd + " > " + output_prefix + "." + suffix)
        else:  # bwa aln
            # alignment step
            log.info("run bwa aln algorithm on %s"%input)
            cmd = "bwa aln {} {} -t {}".format(self.index_prefix, input, thread)
            if max_diff is not None:
                cmd += " -n {} ".format(max_diff)
            if params:
                cmd += " " + params
            cmd += " > {}.sai".format(output_prefix)
            # run alignment step
            self._check_call(cmd)

            # samse step
            cmd = "bwa samse {} {}.sai {}".format(self.index_prefix, output_prefix, input)
            if samse_params:
                cmd += " " + samse_params
            if bam:
                cmd += " | samtools view -bh"
            # run samse step
            suffix = "sam" if not bam else "bam"
            cmd = cmd + " > " + output_prefix + "." + suffix
            self._check_call(cmd)
            os.remove(output_prefix+".sai")
    
    def __del__(self):
        """ close the log file. """
        if self.log_file is not sys.stdout:
            self.log_file.close()
