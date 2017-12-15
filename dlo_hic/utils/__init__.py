from .reverse_complement import reverse_complement

def read_args(args, scope):
    for var in dir(args):
        if not var.startswith("_"):
            scope[var] = getattr(args, var)