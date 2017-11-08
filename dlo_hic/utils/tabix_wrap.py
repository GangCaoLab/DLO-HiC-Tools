import sys
import subprocess

from .parse_text import parse_line_bed6

def sort_bed6(file_in, file_out):
    """ sort bed file. """
    cmd = "sort -k1,1V -k2,2n -k3,3n -u {} > {}".format(file_in, file_out)
    subprocess.check_call(cmd, shell=True)

sort_bedpe_reads1 = sort_bed6

def index_bed6(bedfile):
    """ build tabix index for bedfile.
    please ensure bed file is sorted """
    cmd = "cat {} | bgzip > {}".format(bedfile, bedfile+'.gz')
    subprocess.check_call(cmd, shell=True)
    cmd = "tabix -s 1 -b 2 -e 3 {}".format(bedfile+'.gz')
    subprocess.check_call(cmd, shell=True)

def query_bed6(bedfile, chr_, start, end):
    """ query to indexed bed file """
    start, end = str(start), str(end)
    cmd = "tabix {} {}:{}-{}".format(bedfile+'.gz', chr_, start, end)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    for line in process.stdout:
        yield parse_line_bed6(line)