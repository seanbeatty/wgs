import os

import pypeliner

from wgs.utils import helpers
from wgs.utils import pdfutils


scripts = os.path.join(
    os.path.realpath(os.path.dirname(__file__)),
    'scripts'
)


def hmmcopy_readcounter(input_bam, output_wig, config):
    chromosomes = ['--chromosomes'] + config['hmmcopy_readcounter']['chromosomes']

    cmd = [
              'python',
              os.path.join(scripts, 'read_counter.py'),
              input_bam,
              output_wig,
              '-w',
              str(config['hmmcopy_readcounter']['w']),
              '-m',
              str(config['hmmcopy_readcounter']['m']),
          ] + chromosomes

    helpers.run_cmd(cmd, output=output_wig)


def calc_corr(input_wig, output_file, output_obj, config):
    cmd = [
        'R',
        '--no-save',
        '--args',
        input_wig,
        config['calc_corr']['gc'],
        config['calc_corr']['map'],
        str(config['calc_corr']['mapcutoff']),
        output_file,
        output_obj,
        '<',
        os.path.join(scripts, 'correctReadCount.R'),
    ]

    pypeliner.commandline.execute(*cmd)


def run_hmmcopy(
        tumour_copy,
        tumour_table,
        output_obj,
        output_segments,
        tumour_table_out,
        sample_id,
        config):
    args = {key: 'NULL' if value is None else value
            for key, value in config['run_hmmcopy'].items()}

    param_list = (
        'm',
        'mu',
        'kappa',
        'e',
        'S',
        'strength',
        'lambda',
        'nu',
        'eta',
        'gamma',
    )

    params = [args[param] for param in param_list]

    cmd = [
              'R',
              '--no-save',
              '--args',
              tumour_copy,
              tumour_table,
              args['normal_copy'],
              args['normal_table'],
              output_segments,
              output_obj,
              sample_id,
              args['normal_table_out'],
              tumour_table_out,
          ] + params + ['<', os.path.join(scripts, 'hmmcopy.R')]

    pypeliner.commandline.execute(*cmd)


def plot_hmm(
        tumour_copy,
        hmmcopy_res,
        correction_plots_dir,
        hmmcopy_plots_dir,
        bias_pdf,
        correction_pdf,
        hmmcopy_pdf):
    helpers.makedirs(correction_plots_dir)
    helpers.makedirs(hmmcopy_plots_dir)

    cmd = [
        'Rscript',
        os.path.join(scripts, 'plot_hmmcopy.R'),
        tumour_copy,
        hmmcopy_res,
        correction_plots_dir,
        bias_pdf,
        hmmcopy_plots_dir,
    ]

    pypeliner.commandline.execute(*cmd)

    correction_pdfs = [os.path.join(correction_plots_dir, f)
                       for f in os.listdir(correction_plots_dir) if f.endswith('.pdf')]
    hmmcopy_pdfs = [os.path.join(hmmcopy_plots_dir, f)
                    for f in os.listdir(hmmcopy_plots_dir) if f.endswith('.pdf')]

    pdfutils.merge_pngs(correction_pdfs, correction_pdf)
    pdfutils.merge_pngs(hmmcopy_pdfs, hmmcopy_pdf)


def annot_hmm(input_segments, output_file, config):
    cmd = [
        'python',
        os.path.join(scripts, 'pygene_annotation.py'),
        '--gene_sets_gtf',
        config['annot_hmm']['gtf'],
        '--outfile',
        output_file,
        '--infile',
        input_segments,
    ]

    pypeliner.commandline.execute(*cmd)
