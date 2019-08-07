import logging
log = logging.getLogger(__name__)


class ServerConfigError(Exception):
    pass


def file_exist(path):
    from os.path import exists, isfile
    return exists(path) and isfile(path)


def bwa_index_exist(index_prefix):
    from os.path import exists
    suffix = ["amb", "bwt", "sa", "ann", "pac"]
    files = [index_prefix + "." + s for s in suffix]
    exists = [exists(f) for f in files]
    all_exists = all(exists)
    return all_exists


def is_supported_genome(gname):
    from os.path import abspath, dirname, join, exists
    here = dirname(abspath(__file__))
    built_in_genome_path = join(here, "../templates/chromosome_length")
    return exists(join(built_in_genome_path, gname+'.txt'))


def filter_genomes(genomes):

    def is_valid(g):
        fields = ['name', 'fasta', 'bwa_index_prefix', 'chromosome_file']
        for f in fields:
            if f not in g:  # check fields exists
                log.warning("Reference genome config lack required field: '{}'".format(f))
                return False
            if f == 'bwa_index_prefix':
                if not bwa_index_exist(g['bwa_index_prefix']):
                    log.warning("BWA aligner index prefix is not valid.")
                    return False
            if (f == 'chromosome_file') and is_supported_genome(g[f]):
                return True
            if f in ['fasta']:
                if not file_exist(g[f]):
                    log.warning("field: " + f + " path: '" + path + "' is not exist.")
                    return False
        return True

    valid_genomes = []
    for g in genomes:
        if is_valid(g):
            valid_genomes.append(g)
        else:
            log.warning("Reference genome is not valid: {}".format(str(g)))
    return valid_genomes


def check_required_fields(conf):
    fields = [
        ('host', ''),
        ('port', ''),
        ('password_hash', "Please specify the password_hash code, "
                          "hash code can be generate using 'dlohicwui hash' command."),
    ]
    for f, info in fields:
        if f not in conf:
            if info:
                log.error(info)
                raise ServerConfigError()
            else:
                log.error("{} field is required.".format(f))
                raise ServerConfigError()


def check_file_exist(conf):
    fields = ['juicer_tools_jar']
    for f in fields:
        if not file_exist(conf[f]):
            log.warning("field: " + f + " path: '" + path + "' is not exist.")


def check_java_env(conf):
    from subprocess import PIPE, Popen
    try:
        p = Popen(['java', '-version'], stdout=PIPE, stderr=PIPE)
        _, _ = p.communicate()
        conf['java_exist'] = True
        p.kill()
    except FileNotFoundError as e:
        log.warning(str(e))
        log.warning('Java environment does not exist. Not support convert to .hic file.')
        conf['java_exist'] = False
        return

    jar = conf['juicer_tools_jar']
    if not jar: return

    p = Popen(['java', '-jar', jar], stdout=PIPE, stderr=PIPE)
    _, _ = p.communicate()
    if p.returncode != 0:
        log.warning("Jar file '{}' is not runable.".format(jar))
        conf['juicer_tools_jar'] = ""
    p.kill()


def check_conf(conf):
    check_required_fields(conf)
    check_file_exist(conf)
    check_java_env(conf)


def fill_optional(conf):
    fields = [
        ('reference_genomes', []),
        ('juicer_tools_jar', "")
    ]
    for f, v in fields:
        if f not in conf:
            conf[f] = v


def parse(path):
    import json
    with open(path) as f:
        conf = json.load(f)
    conf['path'] = path
    conf['reference_genomes'] = filter_genomes(conf['reference_genomes'])
    import os
    from os.path import abspath
    conf['current_dir'] = abspath(os.curdir)

    fill_optional(conf)
    check_conf(conf)

    return conf

