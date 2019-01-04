import os
from os.path import join, split, splitext, dirname
import logging
from collections import OrderedDict

import click
from mako.template import Template
from mako.lookup import TemplateLookup

from dlo_hic.utils.pipeline import qc_files, DIRS


log = logging.getLogger(__name__)


def get_sample_ids(pipe_workdir):
    # suppose the pairs.gz file which in the subdir 4 exist
    pairs_step_id = 4
    guess_exist_type = (pairs_step_id, 'pairs.gz')

    dir_ = join(pipe_workdir, DIRS[guess_exist_type[0]])
    print(dir_)
    files = os.listdir(dir_)

    def extract_id(path):
        fname = split(path)[1]
        id_ = fname.replace('.'+guess_exist_type[1], '')
        return id_

    files_ = filter(lambda f: f.endswith(guess_exist_type[1]), files)
    sample_ids = list(set( map(extract_id, files_) ))
    return sample_ids


def get_reads_counts_qc(pipe_workdir, sample_id):
    res = OrderedDict()
    reads_counts_qc = []
    for step, val in qc_files(sample_id).items():
        if step == 'noise_reduce':
            reads_counts_qc.append( (step, val['normal']) )
        elif step == 'bedpe2pairs':
            reads_counts_qc.append( (step, val['counts']) )
        else:
            reads_counts_qc.append( (step, val) )
    for step, qc_file in reads_counts_qc:
        qc_path = join(pipe_workdir, qc_file)
        res[step] = load_reads_counts_qc(qc_path)
    return res


def load_reads_counts_qc(path):
    res = OrderedDict()
    with open(path) as f:
        for line in f:
            line = line.strip()
            items = line.split()
            res[items[0]] = items[1]
    return res


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


def get_qc_contents(pipe_workdir, sample_id):
    reads_counts_qc = get_reads_counts_qc(pipe_workdir, sample_id)

    csv_path = qc_files(sample_id)['bedpe2pairs']['chr_interactions']
    csv_path = join(pipe_workdir, csv_path)
    chr_interactions = load_chr_interaction_csv(csv_path)

    res = {
        'reads_counts': reads_counts_qc,
        'chr_interactions': chr_interactions,
    }
    return res


def render_report(sample_id, qc_contents, report_format):
    if report_format == 'txt':
        return render_txt_report(sample_id, qc_contents)
    else:
        return render_html_report(sample_id, qc_contents)


def render_txt_report(sample_id, qc_contents):
    res = ""
    for step, qc in qc_contents.items():
        res += "[{}]\n".format(step)
        for item, val in qc.items():
            res += "{}\t{}\n".format(item, val)
        res += "\n"
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
    report = template.render(sample_id=sample_id, qc_contents=qc_contents)
    return report


@click.command(name="gen_qc_report")
@click.argument("pipe-workdir")
@click.argument("output")
@click.option("--out-format",
    default="html",
    type=click.Choice(['html', 'txt']),
    help="The format of quility control report, 'html' or 'txt'")
def _main(pipe_workdir, output, out_format):
    for s_id in get_sample_ids(pipe_workdir):
        qc_contents = get_qc_contents(pipe_workdir, s_id)
        log.info("Generating {} format quility control report of sample '{}'.".format(out_format, s_id))
        report = render_report(s_id, qc_contents, out_format)
        with open(output, 'w') as f:
            f.write(report)
            log.info("Quility control report of sample '{}' generated, saving to {}".format(s_id, output))
        print(output)


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
