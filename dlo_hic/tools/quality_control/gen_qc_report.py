import os
from os.path import join, split, splitext, dirname
import logging
from collections import OrderedDict

import click
from mako.template import Template
from mako.lookup import TemplateLookup

import dlo_hic
from dlo_hic.utils.pipeline import qc_files, DIRS


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


def get_qc_contents(pipe_workdir, sample_id):
    load_funcs = OrderedDict({
        'extract_PET': {
            'main': load_pet_main,
            'adapter': load_comp,
            'adapter_svg': load_svg,
        },
        'build_bedpe': {
            'main': load_bedpe_main,
        },
        'noise_reduce': {
            'main': load_comp,
        },
        'bedpe2pairs' : {
            'comp': load_comp,
            'chr_interactions': load_chr_interaction_csv,
            'pet_span_stats': load_pet_span_stats,
            'pet_span_fig': load_svg,
        }
    })
    res = OrderedDict()
    files = qc_files(sample_id)

    qc_dir = join(pipe_workdir, DIRS['qc'])
    if sample_id + '.adapter.txt' in os.listdir(qc_dir):  # if adapter file in path
        files['extract_PET'].update({
            'adapter': join(qc_dir, sample_id + '.adapter.txt'),
            'adapter_svg': join(qc_dir, sample_id + '.adapter.svg'),
        })

    for step in load_funcs:
        res[step] = OrderedDict()
        for item in load_funcs[step]:
            func = load_funcs[step][item]
            file_ = os.path.join(pipe_workdir, files[step][item])
            res[step][item] = func(file_)
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
#        import ipdb; ipdb.set_trace()
        report = render_html_report(s_id, qc_contents)
        with open(output, 'w') as f:
            f.write(report)
            log.info("Quility control report of sample '{}' generated, saving to {}".format(s_id, output))


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
