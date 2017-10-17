from os.path import abspath, join, dirname

here = abspath(__file__)
here_dir = dirname(here)

suffix = '.txt'

supported_genomes = {
    'hg38': join(here_dir, 'hg38' + suffix),
    'hg19': join(here_dir, 'hg19' + suffix),
    'hg18': join(here_dir, 'hg18' + suffix),
    'mm10': join(here_dir, 'mm10' + suffix),
    'mm9' : join(here_dir, 'mm9'  + suffix),
    'mm8' : join(here_dir, 'mm8'  + suffix),
}
