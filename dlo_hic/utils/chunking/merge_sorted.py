from subprocess import Popen, PIPE, check_call


def merge_sorted_pairs(sorted_files, output, ncpu=8):
    order_args = "-k2,2 -k4,4 -k3,3n -k5,5n -k6,6 -k7,7".split()
    cmd = ['sort', '-m'] + order_args + sorted_files + ['--parallel={}'.format(ncpu)]
    cmd = " ".join(cmd) + " > " + output
    check_call(cmd, shell=True)
    return output
