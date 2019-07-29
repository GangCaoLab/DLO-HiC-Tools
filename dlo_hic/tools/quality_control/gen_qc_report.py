import re
import os
from os.path import join, split, splitext, dirname
import logging
from collections import OrderedDict

import click
from mako.template import Template
from mako.lookup import TemplateLookup

import dlo_hic
from dlo_hic.utils.pipeline import qc_files, DIRS, parse_config


log = logging.getLogger(__name__)

pairs_step_id = 4


def get_sample_ids(pipe_workdir):
    """
    Get the sample ids in a pipeline work dir.
    """
    # suppose the pairs.gz file which in the subdir 4 exist
    guess_exist_type = (pairs_step_id, 'pairs.gz')

    dir_ = join(pipe_workdir, DIRS[guess_exist_type[0]])
    files = os.listdir(dir_)

    def extract_id(path):
        fname = split(path)[1]
        id_ = fname.replace('.'+guess_exist_type[1], '')
        return id_

    files_ = filter(lambda f: f.endswith(guess_exist_type[1]), files)
    sample_ids = list(set( map(extract_id, files_) ))
    return sample_ids



def load_chr_interaction_csv(path):
    import pandas as pd
    import numpy as np
    df = pd.read_csv(path, index_col=0)

    data_objs = []
    for idx, value in np.ndenumerate(df.values):
        i, j = idx
        if j >= i:
            obj = {
                'pos': [i, j],
                'chr': [df.index[i], df.index[j]],
                'value': value
            }
            data_objs.append(obj)

    res = {
        'chromosomes': list(df.columns),
        'data': data_objs,
    }
    return res


def load_pet_span_stats(path):
    stats = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("[Describe]"):
                continue
            elif len(line) == 0 or line.startswith("["):
                break
            else:
                items = line.split("\t")
                key = items[0]
                val = float(items[1])
                stats.append((key, val)) 
    return stats


def load_svg(path):
    with open(path) as f:
        contents = f.read()
    contents = contents.strip()
    return contents


def read_tab_splited(fh):
    """
    read a tab splited section in a file, return when read line is empty
    """
    res = OrderedDict()
    for line in fh:
        line = line.strip()
        if line.startswith("#"):
            continue
        elif line == '':
            break
        else:
            items = line.split("\t")
            res[items[0]] = items[1]
    return res


def load_comp(path):
    with open(path) as f:
        content = read_tab_splited(f)
    return content


def load_pet_main(path):
    res = OrderedDict()
    with open(path) as f:
        res['linkers'] = read_tab_splited(f)
        res['flag_stat'] = read_tab_splited(f)
        res['PET1_len_dist'] = read_tab_splited(f)
        res['PET2_len_dist'] = read_tab_splited(f)
        res['linker_match_score_dist'] = read_tab_splited(f)
        res['adapter_match_score_dist'] = read_tab_splited(f)
    return res


def load_bedpe_main(path):
    res = OrderedDict()
    with open(path) as f:
        res['PET1'] = read_tab_splited(f)
        res['PET2'] = read_tab_splited(f)
        res['paired'] = read_tab_splited(f)
    return res


def load_noise_reduce_main(path):
    res = OrderedDict()
    with open(path) as f:
        res['type'] = read_tab_splited(f)
        res['position'] = read_tab_splited(f)
        res['distance_PET1'] = read_tab_splited(f)
        res['distance_PET2'] = read_tab_splited(f)
    return res


def load_timepoints(log_file):
    """ Extract key points timestamp from log file. """
    from datetime import datetime
    from dlo_hic.config import LOGGING_DATE_FMT
    key_time_points = OrderedDict()
    with open(log_file) as f:
        for line in f:
            line = line.strip()
            match = re.search("@ (.*): (START|END): (.*)", line)
            if match:
                time_, flag, pname = match.groups()
                time = datetime.strptime(time_, LOGGING_DATE_FMT)
                key_time_points.setdefault(pname, dict())
                key_time_points[pname][flag.lower()] = time
    for _, times in key_time_points.items():
        start = times.get('start')
        end = times.get('end')
        diff = end - start if start and end else None
        times['start'] = str(start) if start else ''
        times['end']   = str(end) if end else ''
        times['diff']  = str(diff) if diff else ''
    return key_time_points


def load_pipe_config(config_file):
    lines = []
    with open(config_file) as f:
        for line in f:
            line = line.strip()
            if re.match("\s*#", line):  # skip comments
                continue
            line = re.sub("#.*$", "", line)
            if line.startswith('['):
                lines.append('')
            if line != '':
                lines.append(line)
    res = "\n".join(lines)
    return res


def get_output_paths(pipe_workdir):
    pass


def get_software_version_info():
    from subprocess import Popen, PIPE
    import re
    s2v = OrderedDict()

    def get_version(command, pattern, *args):
        try:
            p = Popen([command] + list(args), stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            return 'Not found.'
        out, err = p.communicate()
        err = err.decode('utf-8')
        out = out.decode('utf-8')
        for line in (err+out).split("\n"):
            m = re.search(pattern, line)
            if m:
                return line.strip()

    for cmd in ['bwa', 'samtools', 'tabix', 'pairix']:
        s2v[cmd] = get_version(cmd, "[Vv]ersion")
    for cmd in ['cooler', 'mafft']:
        s2v[cmd] = get_version(cmd, ".*", "--version")

    res = "\n".join([s+": "+v for s,v in s2v.items()])
    
    return res


def get_python_package_info():
    from subprocess import Popen, PIPE
    p = Popen(['pip', 'freeze'], stdout=PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    return out


def split_iter(iterator, f):
    yes_list = []
    no_list = []
    for i in iterator:
        if f(i):
            yes_list.append(i)
        else:
            no_list.append(i)
    return yes_list, no_list


def get_qc_contents(pipe_workdir, sample_id):
    """
    Compose a dict for render qc report.
    """

    load_funcs = OrderedDict({
        'extract_PET': {
            'main': load_pet_main,
        },
        'build_bedpe': {
            'main': load_bedpe_main,
        },
        'noise_reduce': {
            'main': load_noise_reduce_main,
        },
        'bedpe2pairs' : {
            'comp': load_comp,
            'chr_interactions': load_chr_interaction_csv,
            'pet_span_stats': load_pet_span_stats,
            'pet_span_fig': load_svg,
        },
        'other' : {
            'time_points': load_timepoints,
            'pipe_config': load_pipe_config,
        }
    })
    res = OrderedDict()
    files = qc_files(sample_id)

    # if infer_adapter has output, update load_funcs and files
    qc_dir = DIRS['qc']
    contents_qc_dir = os.listdir(join(pipe_workdir, qc_dir))
    if sample_id + '.adapter.txt' in contents_qc_dir:
        load_funcs['extract_PET'].update({
            'adapter': load_comp,
            'adapter_svg': load_svg,
        })
        files['extract_PET'].update({
            'adapter': join(qc_dir, sample_id + '.adapter.txt'),
            'adapter_svg': join(qc_dir, sample_id + '.adapter.svg'),
        })

    # if global cmap is ploted, add it to files
    if sample_id + '.global_cmap.svg' in contents_qc_dir:

        def load_chroms_svg(path):
            contents = []
            for p in os.listdir(path):
                if not p.endswith('.svg'):
                    continue
                chrom = splitext(p)[0]
                svg = load_svg(join(path, p))
                contents.append((chrom, svg))

            num_chrs, nonnum_chrs = split_iter(contents, lambda t: t[0].replace('chr', '').isdigit())
            contents = sorted(num_chrs, key=lambda t:int(t[0].replace('chr', ''))) + sorted(nonnum_chrs)
            chroms2svg = OrderedDict(contents)
            return chroms2svg

        load_funcs.update({
            'contact_map' : {
                'global': load_svg,
                'chroms': load_chroms_svg  # load all imgs in path
            }
        })
        files.update({
            'contact_map': {
                'global': join(qc_dir, sample_id + '.global_cmap.svg'),
                'chroms': join(qc_dir, sample_id + '.chroms_cmap')
            }
        })

    # add log file to files
    log_dir = DIRS['log']
    files['other'] = {'time_points': join(log_dir, sample_id + '.log')}

    # add pip config file to files
    files['other'].update({'pipe_config': join(qc_dir, sample_id + '.pipe_config.ini')})

    # load contents from files
    for step in load_funcs:
        res[step] = OrderedDict()
        for item in load_funcs[step]:
            func = load_funcs[step][item]
            file_ = os.path.join(pipe_workdir, files[step][item])
            res[step][item] = func(file_)

    # load dependency info
    res.update({
        'dependency': {
            'software': get_software_version_info(),
            'python': get_python_package_info(),
        }
    })

    return res


def render_html_report(sample_id, qc_contents):
    here = os.path.dirname(os.path.abspath(__file__))
    template_path = join(here, "../../templates/qc_report/qc_report.mako")
    template_dir = dirname(template_path)
    lookup = TemplateLookup(directories=[template_dir],
                            input_encoding='utf-8',
                            output_encoding='utf-8',
                            encoding_errors='replace')
    template = Template(filename=template_path, lookup=lookup)
    from datetime import datetime
    now = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%m:%d')
    report = template.render(sample_id=sample_id, qc_contents=qc_contents, version=dlo_hic.__version__, datetime=now)
    return report


@click.command(name="gen_qc_report")
@click.argument("pipe-workdir")
@click.argument("sample-id")
@click.argument("output")
@click.option("--out-format",
    default="html",
    type=click.Choice(['html', 'txt']),
    help="The format of quility control report, 'html' or 'txt'")
def _main(pipe_workdir, sample_id, output, out_format):
    """ Generate pipeline quality control report. """
    possible_ids = get_sample_ids(pipe_workdir)
    s_id = sample_id
    if s_id not in possible_ids:
        raise IOError("{} is not a valid sample id, candidate sample ids: {}".format(s_id, possible_ids))

    log.info("Generating {} format quility control report of sample '{}'.".format(out_format, s_id))
    if out_format == 'txt':
        import subprocess as subp
        cmd = "tail -n +1 {}/*.txt > {}".format(os.path.join(pipe_workdir, DIRS['qc']), output)
        subp.check_call(cmd, shell=True)
    else:
        qc_contents = get_qc_contents(pipe_workdir, s_id)
        report = render_html_report(s_id, qc_contents)
        with open(output, 'w') as f:
            f.write(report)
            log.info("Quility control report of sample '{}' generated, saving to {}".format(s_id, output))


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
