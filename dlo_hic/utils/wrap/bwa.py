import os
import sys
from datetime import datetime
import logging
import subprocess
from subprocess import Popen, PIPE
from io import TextIOWrapper
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

    outfmt : {'sam', 'bam'}
        output format.
    
    log_file : {str, '-'}
        path to log file, use '-' for output log to stdout.
    """
    def __init__(self, index_prefix, algorithm='aln', outfmt='bam', log_file='-'):
        self.index_prefix = index_prefix
        self.algorithm = algorithm
        self.log_file = log_file
        self.outfmt = outfmt

    def _check_call(self, cmd):
        """ 
        call cmd and log to self.log_file
        """
        process = Popen(cmd, stderr=PIPE, shell=True)
        if self.log_file != '-':
            log_f = open(self.log_file, 'a')
        with process.stderr:
            for line in process.stderr:
                line = line.decode("utf-8")
                timestr = datetime.now().strftime("%m/%d/%y %H:%M:%S")
                line = "[{}]".format(timestr) + line # add datetime information
                log_f.write(line)
        if self.log_file != '-':
            log_f.close()
        exitcode = process.wait() 
        if exitcode != 0:
            raise subprocess.CalledProcessError(exitcode, cmd)

    def index(self, ref_fasta):
        """ build index """
        index_prefix = self.index_prefix
        msg = "build bwa index. fasta file: %s, index_prefix: %s"%(ref_fasta, index_prefix)
        log.info(msg)
        cmd = "bwa index -p {} {}".format(index_prefix, ref_fasta)
        self._check_call(cmd)

    def run(self, input, output_prefix,
            thread=cores,
            max_diff=None,
            params=None):
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
            if self.outfmt == "bam":
                cmd += " | samtools view -bh"

            # run cmd
            suffix = "sam" if outfmt == "sam" else "bam"
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
            self._output_prefix = output_prefix
            self._input_fq = input
            # run alignment step
            self._check_call(cmd)

    def samse(self, samse_params=None, stdout=False):
        if self.algorithm != 'aln':
            log.warning("Only aln algorithm need samse.")
            return

        # samse step
        cmd = "bwa samse {} {}.sai {}".format(self.index_prefix, self._output_prefix, self._input_fq)
        if samse_params:
            cmd += " " + samse_params
        if stdout:
            if self.log_file == '-':
                if self.outfmt == "bam":
                    cmd += " | samtools view -bh"
                p = Popen(cmd, stdout=PIPE, shell=True)
            else:
                cmd = cmd + " 2>> " + self.log_file
                if self.outfmt == "bam":
                    cmd += " | samtools view -bh"
                p = Popen(cmd, stdout=PIPE, shell=True)
            return TextIOWrapper(p.stdout)
        else:
            if self.outfmt == "bam":
                cmd += " | samtools view -bh"
            suffix = "sam" if not self.outfmt == "sam" else "bam"
            cmd = cmd + " > " + self._output_prefix + "." + suffix
            self._check_call(cmd)
            os.remove(self._output_prefix+".sai")
    