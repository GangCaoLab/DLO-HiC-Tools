"""
utils for pipeline configuration.
"""

import logging
import os
from os.path import join
import re
from collections import OrderedDict
import subprocess

from .parse_text import is_comment


DIRS = {
    "qc": "00-qc",
   "log": "00-log",
       1: "01-extpet",
       2: "02-bedpe",
       3: "03-nr",
       4: "04-rr",
       5: "05-pairs",
       6: "06-result",
}


OUTPUT_FILE_TYPES = [
    ["pet.fq"],
    ["pet.bam", "pet.filtered.bam", "pet.bed", "uniq.bedpe"],
    ["nr.bedpe", "nr.bedpe.err",],
    ["rr.bedpe"],
    ["pairs", "pairs.gz", "pairs.gz.px2"],
    ["hic", "cool"],
]


def make_result_dirs():
    for d in DIRS.values():
        if not os.path.exists(d):
            os.mkdir(d)


class SnakeFilter(logging.Filter):
    """
    Log filter for filter out 'snakemake.logging' logger's output
    """
    def filter(self, record):
        if record.name == "snakemake.logging":
            return False
        else:
            return True


def parse_config(config_file):
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = config2dict(config)
    check_config(config_dict)
    return config_dict


def config2dict(config):
    """
    Parse config dict, convert string fields to correct type.
    """
    parsed_config = {}
    for sec_name, section in config.items():
        parsed_config[sec_name] = {}
        for field, val in section.items():
            if val.strip() and not is_comment(val):
                parsed_val = eval(val)
            else:
                parsed_val = None
            parsed_config[sec_name][field] = parsed_val
    return parsed_config


def check_config(config):
    """
    Check the parsed config.
    """
    check_required(config)
    check_required(config)

    ## other checking
    # check result JuicerToolsJar if ResultFormats contain '.hic'
    if '.hic' in config['RESULT']['resultformats']:
        if not config['RESULT']['juicertoolsjar']:
            raise IOError("RESULT/juicertoolsjar must be specified if resultformats contain '.hic'")


def check_required(config):
    """ check required fields"""
    required_fields = [
        'GLOBAL/loglevel',
        'DATA/input_dir',
        'DATA/fasta',
        'DATA/bwaindexprefix',
        'DATA/restriction_site',
        'DATA/restriction_name',
        'EXTRACT_PET/linker-a',
    ]

    for require in required_fields:
        section, field = require.split("/")
        if not config[section][field]:
            raise IOError("%s/%s is required."%(section, field))


def check_files(config):
    """ check files exist or not. """
    file_fields = [
        'GLOBAL/workingdir',
        'DATA/input_dir',
        'DATA/fasta',
        'NOISE_REDUCE/restriction_sites_bed',
        'RESULT/juicertoolsjar',
    ]

    for file_ in file_fields:
        section, field = file_.split("/")
        if config[section][field]: # if field specified check path exist or not
            if not os.path.exists(config[section][field]):
                raise IOError("Path to %s/%s is not exist."%(section, field))


class PipelineSetting():
    pass


def load_global_setting(config):
    setting = PipelineSetting()
    setting.ncpu = config['GLOBAL']['numbercpus']
    setting.input_dir = config['DATA']['input_dir']
    setting.is_qc = config['PROCESSES']['qc']
    setting.qc_report_prefix = config['QUALITY_CONTROL']['qc_report_prefix']
    setting.qc_report_suffix = config['QUALITY_CONTROL']['report_format']
    setting.working_dir = config['GLOBAL']['workingdir'] or "./"
    setting.result_formats = config['RESULT']['resultformats']
    setting.keep = config['PROCESSES']['keep']
    return setting


def get_input_fastq(input_dir):
    def get_input_fastq_(wildcard):
        possible = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
        for suffix in possible:
            path = os.path.join(input_dir, wildcard.sample + suffix)
            if os.path.exists(path):
                return path
    return get_input_fastq_


def fetch_all_fastq(input_dir):
    result = []
    for f in os.listdir(input_dir):
        if f.endswith(".fq") or \
           f.endswith(".fastq") or \
           f.endswith(".fq.gz") or \
           f.endswith(".fastq.gz"):
           result.append(os.path.join(input_dir, f))
    return result


def get_samples_id(fastq_files):
    result = []
    for fq in fastq_files:
        f = os.path.basename(fq)
        sid = re.sub("\.fastq.gz$|\.fq.gz$|\.fastq$|\.fq$", "", f)
        result.append(sid)
    return result


def local_logger(wildcard):
    from dlo_hic.config import LOGGING_DATE_FMT, LOGGING_FMT
    sample = wildcard.sample
    logger = logging.getLogger("pipeline-"+sample)

    if not logger.handlers:
        log_file = DIRS['log'] + "/" + sample + ".log"
        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter(
            fmt=LOGGING_FMT,
            datefmt=LOGGING_DATE_FMT)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        def log_stage_boundary(message, separator='=', num=30):
            """
            log a stage boundary.
            """
            logger.info(separator*num)
            logger.info(message)
            logger.info(separator*num)
        logger.log_stage_boundary = log_stage_boundary

    return logger


def sub_dir(id_):
    return DIRS[id_]


def qc_logs(sample_id):
    return OrderedDict([
        ('extract_PET',       join(sub_dir(1), sample_id + '.qc.pet.txt')),
        ('build_bedpe',       join(sub_dir(2), sample_id + '.qc.bedpe.txt')),
        ('noise_reduce',      join(sub_dir(3), sample_id + '.qc.nr.txt')),
        ('noise_reduce.err',  join(sub_dir(3), sample_id + '.qc.nr.err.txt')),
        ('remove_redundancy', join(sub_dir(4), sample_id + '.qc.rr.txt')),
    ])


def output_files(setting):
    def is_keep(k):
        if setting.keep == 'ALL':
            return True
        else:
            return k in setting.keep

    def temp_file(fname, key_):
        """
        mark file as temp file, if key_ not in config['PROCESSES']['keep']
        """
        from snakemake.rules import temp
        if is_keep(key_):
            return fname
        else:
            return temp(fname)

    def output_files_(sample_id):

        res = OrderedDict([])

        def path(dir, sample, ext):
            return join(dir, sample + '.' + ext)

        for i in range(len(OUTPUT_FILE_TYPES)):
            dir_ = DIRS[i+1]
            for ext in OUTPUT_FILE_TYPES[i]:
                if 'pet' in ext:
                    ext1 = ext.replace('pet', 'pet1')
                    ext2 = ext.replace('pet', 'pet2')
                    res[ext1] = path(dir_, sample_id, ext1)
                    res[ext2] = path(dir_, sample_id, ext2)
                else:
                    res[ext] = path(dir_, sample_id, ext)
        
        for k_, p_ in res.items():
            res[k_] = temp_file(p_, k_)

        return res
    
    return output_files_


def get_targets(setting):
    from snakemake.rules import expand

    all_fastq = fetch_all_fastq(setting.input_dir)
    all_sample = get_samples_id(all_fastq)

    all_qc_report = expand(join(DIRS["qc"],  setting.qc_report_prefix+"{sample}."+setting.qc_report_suffix), sample=all_sample) if setting.is_qc else []
    all_pairs     = expand(join(DIRS[5], "{sample}.pairs.gz"), sample=all_sample)
    all_hic       = expand(join(DIRS[6], "{sample}.hic"),      sample=all_sample) if '.hic' in setting.result_formats else []
    all_cool      = expand(join(DIRS[6], "{sample}.cool"),     sample=all_sample) if '.cool' in setting.result_formats else []
    all_ = all_pairs + all_hic + all_cool + all_qc_report

    return all_


def chromosome_files():
    here = os.path.dirname(os.path.abspath(__file__))
    template_path = join(here, "../templates/chromosome_length")
    files = os.listdir(template_path)
    file2id = lambda fname: fname.split(".")[0]
    res = {
        file2id(f): join(template_path, f) for f in files
    }
    return res


def gen_chromosome_file(genomeid, out_path):
    genomeid2path =  chromosome_files()
    chromosomes = list(genomeid2path.keys())
    if genomeid not in genomeid2path:
        msg = "If you are not use {}, you must ".format("/".join(chromosomes)) + \
              "specify the path to the chromosome file."
        raise ValueError(msg)
    if not os.path.exists(out_path):
        subprocess.check_call(["cp", genomeid2path[genomeid], out_path])
