"""
utils for pipeline configuration.
"""

from datetime import datetime
import logging
import os

from .parse_text import is_comment


__all__ = ['Timmer', 'SnakeFilter', 'parse_config', 'check_config']


class Timmer(object):
    """
    Context manager for record key time points and delta.
    """
    def __init__(self, name, recorder):
        """
        :name: process name
        :recorder: OrderedDict or filename,
            where time point and delta to store.
        """
        self.name = name
        self.recorder = recorder

    def __enter__(self):
        self.before = datetime.now()
        if isinstance(self.recorder, dict):
            self.recorder['before-'+self.name] = self.before
        else:
            self.record_file = open(self.recorder, 'a')
            oline = "before-{}\t{}\n".format(self.name, str(self.before))
            self.record_file.write(oline)

    def __exit__(self, exc_type, value, trackback):
        self.after = datetime.now()
        self.cost = self.after - self.before
        if isinstance(self.recorder, dict):
            self.recorder['after-'+self.name] = self.after
            self.recorder['cost-'+self.name] = self.cost
        else:
            oline = "after-{}\t{}\n".format(self.name, str(self.after))
            self.record_file.write(oline)
            oline = "cost-{}\t{}\n".format(self.name, str(self.cost))
            self.record_file.write(oline)
            self.record_file.close()


class SnakeFilter(logging.Filter):
    """
    Log filter for filter out 'snakemake.logging' logger's output
    """
    def filter(self, record):
        if record.name == "snakemake.logging":
            return False
        else:
            return True


def parse_config(config):
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
    ## check required fields
    required_fields = [
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

    ## check files exist
    file_fields = [
        'DATA/input_dir',
        'DATA/fasta',
        'GLOBAL/workingdir',
        'BUILD_BEDPE/bwa_log',
        'NOISE_REDUCE/restriction_sites_bed',
        'RESULT/juicertoolsjar',
    ]

    for file_ in file_fields:
        section, field = file_.split("/")
        if config[section][field]: # if field specified check path exist or not
            if not os.path.exists(config[section][field]):
                raise IOError("Path to %s/%s is not exist."%(section, field))

    ## other checking
    # check result JuicerToolsJar if ResultFormats contain '.hic'
    if '.hic' in config['RESULT']['resultformats']:
        if not config['RESULT']['juicertoolsjar']:
            raise IOError("RESULT/juicertoolsjar must be specified if resultformats contain '.hic'")
