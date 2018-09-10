import os
from os.path import join, split, splitext, dirname
import logging
from collections import OrderedDict

import click
from mako.template import Template
from mako.lookup import TemplateLookup

from dlo_hic.utils.pipeline import qc_logs, sub_dir


log = logging.getLogger(__name__)


def get_sample_ids(pipe_workdir):
    # suppose the pairs.gz file which in the subdir 5 exist
    guess_exist_type = (5, 'pairs.gz')

    dir_ = join(pipe_workdir, sub_dir(guess_exist_type[0]))
    files = os.listdir(dir_)

    def extract_id(path):
        fname = split(path)[1]
        id_ = fname.replace('.'+guess_exist_type[1], '')
        return id_

    files_ = filter(lambda f: f.endswith(guess_exist_type[1]), files)
    sample_ids = list(set( map(extract_id, files_) ))
    return sample_ids


def get_qc_contents(pipe_workdir, sample_id):
    res = OrderedDict()
    for step, qc_file in qc_logs(sample_id).items():
        qc_path = join(pipe_workdir, qc_file)
        res[step] = load_qc(qc_path)
    return res


def load_qc(path):
    res = OrderedDict()
    with open(path) as f:
        for line in f:
            line = line.strip()
            items = line.split()
            res[items[0]] = items[1]
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
@click.argument("output-dir")
@click.option("--out-format",
    default="html",
    help="The format of quility control report, 'html' or 'txt'")
def _main(pipe_workdir, output_dir, out_format):
    for s_id in get_sample_ids(pipe_workdir):
        qc_contents = get_qc_contents(pipe_workdir, s_id)
        log.info("Generating {} format quility control report of sample '{}'.".format(out_format, s_id))
        report = render_report(s_id, qc_contents, out_format)
        output = join(output_dir, s_id + "." + out_format)
        with open(output, 'w') as f:
            f.write(report)
            log.info("Quility control report of sample '{}' generated, saving to {}".format(s_id, output))
        print(output)


main = _main.callback


if "__name__" == "__main__":
    _main()
