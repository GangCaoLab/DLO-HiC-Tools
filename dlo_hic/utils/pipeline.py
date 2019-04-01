"""
utils for pipeline configuration.
"""

import logging
import os
from os.path import join, abspath
import re
from collections import OrderedDict
import subprocess
from functools import partial

from .parse_text import is_comment


log = logging.getLogger(__name__)


# sub directories(under pipeline working dir), for store results
DIRS = {
    "qc": "0a-qc",
   "log": "0b-log",
       1: "01-extpet",
       2: "02-bedpe",
       3: "03-noise_reduced",
       4: "04-valid",
       5: "05-matrix",
}

OUTPUT_FILE_TYPES = [
    ["pet.fq", "trim_flag"],
    ["pet.bam", "pet.uniq.bam", "pet.mul.bam", "pet.unm.bam", "pet.bed", "uniq.bedpe"],
    ["nr.bedpe", "nr.bedpe.sel", "nr.bedpe.re"],
    ["pairs", "pairs.gz", "pairs.gz.px2"],
    ["hic", "cool", "mcool"],
]

# These files store file path
FILE_FIELDS = [
    'GLOBAL/working_dir',
    'DATA/input_dir',
    'DATA/fasta',
    'NOISE_REDUCE/restriction_sites_file',
    'RESULT/juicer_tools_jar',
]


def make_result_dirs(workdir=None):
    """
    Make result directories, like 01-extpet/ 02-bedpe/ ...
    """
    if workdir:
        dirs = {k: join(workdir, dir_) for k, dir_ in DIRS.items()}
    else:
        dirs = DIRS
    for d in dirs.values():
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
    """
    Parse config file,
    return a dict which store pipeline configurations.
    """
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = config2dict(config)
    to_abs_path(config_dict)
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
            val = val.strip()
            if val and not is_comment(val):
                try:
                    parsed_val = eval(val)
                except Exception as e:
                    log.error("Detected some error when parsing the pipeline config file.")
                    log.error("Error occured when parsing field: {}".format(field))
                    raise e
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
    if '.hic' in config['RESULT']['result_formats']:
        if not config['RESULT']['juicer_tools_jar']:
            raise IOError("RESULT/juicertoolsjar must be specified if resultformats contain '.hic'")


def check_required(config):
    """ check required fields """
    required_fields = [
        'GLOBAL/log_level',
        'DATA/input_dir',
        'DATA/fasta',
        'DATA/bwa_index_prefix',
        'DATA/restriction_site',
        'DATA/restriction_name',
        'DATA/chromosome_file',
        'EXTRACT_PET/linker_a',
    ]

    for require in required_fields:
        section, field = require.split("/")
        if not config[section][field]:
            raise IOError("%s/%s is required."%(section, field))


def to_abs_path(config):
    """
    convert paths in the dict to absolute path.
    """
    for file_ in FILE_FIELDS:
        section, field = file_.split("/")
        if config[section][field]:
            # if field specified convert to abs path
            config[section][field] = abspath(config[section][field])


def check_files(config):
    """ check files exist or not. """
    for file_ in FILE_FIELDS:
        section, field = file_.split("/")
        if config[section][field]: # if field specified check path exist or not
            if not os.path.exists(config[section][field]):
                raise IOError("Path to %s/%s is not exist."%(section, field))


class PipelineSetting():
    pass


def load_global_setting(config):
    setting = PipelineSetting()
    setting.ncpu = config['GLOBAL']['number_cpus']
    setting.input_dir = config['DATA']['input_dir']
    setting.chromosome_file = config['DATA']['chromosome_file']
    setting.is_qc = config['PROCESSES']['is_qc']
    setting.qc_report_format = config['QUALITY_CONTROL']['report_format']
    setting.working_dir = config['GLOBAL']['working_dir'] or "./"
    setting.result_formats = config['RESULT']['result_formats']
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


def local_logger(setting):
    def local_logger_(wildcard):
        from dlo_hic.config import LOGGING_DATE_FMT, LOGGING_FMT
        sample = wildcard.sample
        logger = logging.getLogger("pipeline-"+sample)

        if not logger.handlers:
            log_file = join(setting.working_dir, DIRS['log'], sample + ".log")
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
    return local_logger_


def sub_dir(setting):
    workdir = setting.working_dir
    def sub_dir_(id_):
        return join(workdir, DIRS[id_])
    return sub_dir_


def qc_files(sample_id, dirs=DIRS):
    """
    Quality control files about one sample.
    """
    dir_ = dirs['qc']
    return OrderedDict([
        ('extract_PET',       {
            'main':                join(dir_, sample_id + '.pet.txt'),
        }),
        ('build_bedpe',       {
            'main':                join(dir_, sample_id + '.bedpe.txt'),
        }),
        ('noise_reduce',      {
            'main':                join(dir_, sample_id + '.nr.txt'),
        }),
        ('bedpe2pairs',       {
            'comp':                join(dir_, sample_id + '.valid.comp.txt'),
            'chr_interactions':    join(dir_, sample_id + '.valid.chr_interactions.csv'),
            'pet_span_stats':      join(dir_, sample_id + '.pet_span.txt'),
            'pet_span_fig':        join(dir_, sample_id + '.pet_span.svg'),
        }),
    ])


def pipeline_qc_files(setting):
    """
    Quality control files path related to pipeline working dir.
    """
    workdir = setting.working_dir
    dirs = {
        k: join(workdir, dir_)
            for k, dir_ in DIRS.items()
    }
    return partial(qc_files, dirs=dirs)


def all_qc_files(sample_id, workdir=None):
    """
    Get all quality files path
    """
    res = []
    for item in qc_files(sample_id).values():
        if isinstance(item, dict):
            for i in item.values():
                path = join(workdir, i) if workdir else i
                res.append(path)
        else:
            path = join(workdir, item) if workdir else item
            res.append(path)
    return res


def output_files(setting):
    """
    Output files for pipeline.
    """
    def is_keep(k):
        if setting.keep == 'ALL':
            return True
        else:
            if k in OUTPUT_FILE_TYPES[-2] or k in OUTPUT_FILE_TYPES[-1]:
                # keep the result files
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
            dir_ = join(setting.working_dir, DIRS[i+1])
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
    output_files_ = output_files(setting)

    all_fastq = fetch_all_fastq(setting.input_dir)
    all_sample = get_samples_id(all_fastq)

    all_qc_report = expand(join(setting.working_dir, DIRS["qc"], "{sample}."+setting.qc_report_format), sample=all_sample) if setting.is_qc else []
    all_pairs     = expand(output_files_("{sample}")['pairs.gz'], sample=all_sample)
    all_hic       = expand(output_files_("{sample}")['hic'],      sample=all_sample) if '.hic' in setting.result_formats else []
    all_cool      = expand(output_files_("{sample}")['cool'],     sample=all_sample) if '.cool' in setting.result_formats else []
    all_ = all_pairs + all_hic + all_cool + all_qc_report

    return all_


def supported_chromosomes():
    here = os.path.dirname(os.path.abspath(__file__))
    template_path = join(here, "../templates/chromosome_length")
    files = os.listdir(template_path)
    file2id = lambda fname: fname.split(".")[0]
    res = {
        file2id(f): join(template_path, f) for f in files
    }
    return res


def copy_chromosome_file(genomeid, out_path):
    """ copy a chromosome file from template """
    genomeid2path =  supported_chromosomes()
    chromosomes = list(genomeid2path.keys())
    if genomeid not in genomeid2path:
        msg = "If you are not use {}, you must ".format("/".join(chromosomes)) + \
              "specify the path to the chromosome file."
        raise ValueError(msg)
    if not os.path.exists(out_path):
        subprocess.check_call(["cp", genomeid2path[genomeid], out_path])
