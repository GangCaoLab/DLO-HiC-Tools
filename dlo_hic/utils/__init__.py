from .reverse_complement import reverse_complement
from .bwa_wrap import BWA

def read_args(args, scope):
    for var in dir(args):
        if not var.startswith("_"):
            scope[var] = getattr(args, var)