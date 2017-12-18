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
    def __init__(self, index_prefix, log_file='-'):
        self.index_prefix = index_prefix
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

    def index(self, ref_fasta):
        """ build index """
        index_prefix = self.index_prefix
        msg = "build bwa index. fasta file: %s, index_prefix: %s"%(ref_fasta, index_prefix)
        log.info(msg)
        cmd = "bwa index -p {} {}".format(index_prefix, ref_fasta)
        self._check_call(cmd)

    def run(self, input, output_prefix, thread=cores, mem=False, bam=True):
        """ launch bwa aln/mem command, return process """
        # check index exist or not
        index_prefix = self.index_prefix
        suffix = ["amb", "bwt", "sa", "ann", "pac"]
        files = [index_prefix + "." + s for s in suffix]
        exists = [os.path.exists(f) for f in files]
        all_exists = all(exists)
        if not all_exists:
            raise IOError('bwa index not found, run BWA.index firstly.')

        if mem:
            log.info("run bwa mem algorithm on %s"%input)
            if not bam:
                aln_cmd = "bwa mem -t {} {} {} > {}.sam".format(
                    thread, self.index_prefix, input, output_prefix)
            else:
                aln_cmd = "bwa mem -t {} {} {} | samtools view -bh > {}.bam".format(
                    thread, self.index_prefix, input, output_prefix)
            self._check_call(aln_cmd)
        else:
            log.info("run bwa backtrack algorithm on %s"%input)
            aln_cmd = "bwa aln -t {} {} {} > {}.sai".format(
                thread, self.index_prefix, input, output_prefix)
            self._check_call(aln_cmd)
            if not bam:
                trans_cmd = "bwa samse {} {}.sai {} > {}.sam".format(
                    self.index_prefix, output_prefix, input, output_prefix)
            else:
                trans_cmd = "bwa samse {} {}.sai {} | samtools view -bh > {}.bam".format(
                    self.index_prefix, output_prefix, input, output_prefix)
            self._check_call(trans_cmd)
            os.remove(output_prefix+'.sai')
    
    def __del__(self):
        if self.log_file is not sys.stdout:
            self.log_file.close()
