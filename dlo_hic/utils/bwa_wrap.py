import os
import subprocess
import multiprocessing

cores = multiprocessing.cpu_count()

class BWA():
    @classmethod
    def index(cls, ref_fasta, index_prefix):
        suffix = ["amb", "bwt", "sa", "ann", "pac"]
        files = [index_prefix + "." + s for s in suffix]
        exists = [os.path.exists(f) for f in files]
        all_exists = all(exists)
        if not all_exists:
            cmd = "bwa index -p {} {}".format(index_prefix, ref_fasta)
            subprocess.check_call(cmd, shell=True)

    def __init__(self, index_prefix):
        self.index_prefix = index_prefix

    def run(self, input, output_prefix, thread=cores, mem=False, bam=True):
        """ launch bwa aln/mem command, return process """
        if mem:
            if not bam:
                aln_cmd = "bwa mem -t {} {} {} > {}.sam".format(
                    thread, self.index_prefix, input, output_prefix)
            else:
                aln_cmd = "bwa mem -t {} {} {} | samtools view -bh > {}.bam".format(
                    thread, self.index_prefix, input, output_prefix)
            subprocess.check_call(aln_cmd, shell=True)
        else:
            aln_cmd = "bwa aln -t {} {} {} > {}.sai".format(
                thread, self.index_prefix, input, output_prefix)
            subprocess.check_call(aln_cmd, shell=True)
            if not bam:
                trans_cmd = "bwa samse {} {}.sai {} > {}.sam".format(
                    self.index_prefix, output_prefix, input, output_prefix)
            else:
                trans_cmd = "bwa samse {} {}.sai {} | samtools view -bh > {}.bam".format(
                    self.index_prefix, output_prefix, input, output_prefix)
            subprocess.check_call(trans_cmd, shell=True)
            os.remove(output_prefix+'.sai')