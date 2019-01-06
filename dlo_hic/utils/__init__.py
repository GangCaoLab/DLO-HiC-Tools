from .reverse_complement import reverse_complement
from .guess_fq_phred import guess_fq_phred

def read_args(args, scope):
    for var in dir(args):
        if not var.startswith("_"):
            scope[var] = getattr(args, var)
