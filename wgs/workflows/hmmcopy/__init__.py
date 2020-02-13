import os

import pypeliner
import pypeliner.managed as mgd
import tasks


def create_hmmcopy_workflow(bam_file, out_dir, config, sample_id):
    bias_plots_pdf = os.path.join(out_dir, 'plots', '{sample_id}_bias.pdf')
    correction_plots_pdf = os.path.join(out_dir, 'plots', '{sample_id}_correction.pdf')
    hmmcopy_plots_pdf = os.path.join(out_dir, 'plots', '{sample_id}_hmmcopy.pdf')
    tumour_table_out = os.path.join(out_dir, '{}_tumour_correctreads_with_state.txt'.format(sample_id))
    pygene_outfile = os.path.join(out_dir, '{}_hmmcopy.seg.pygenes'.format(sample_id))

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='hmmcopy_readcounter',
        func=tasks.hmmcopy_readcounter,
        args=(
            mgd.InputFile('input.bam', fnames=bam_file),
            mgd.TempOutputFile('infile.wig'),
            config,
        )
    )

    workflow.transform(
        name='calc_corr',
        func=tasks.calc_corr,
        args=(
            mgd.TempInputFile('infile.wig'),
            mgd.TempOutputFile('infile_copy.txt'),
            mgd.TempOutputFile('infile_copy.obj'),
            config,
        )
    )

    workflow.transform(
        name='run_hmmcopy',
        func=tasks.run_hmmcopy,
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('infile_copy.txt'),
            mgd.TempOutputFile('hmmcopy_res.obj'),
            mgd.TempOutputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(tumour_table_out),
            mgd.InputInstance('sample_id'),
            config,
        )
    )

    workflow.transform(
        name='plot_hmm',
        func=tasks.plot_hmm,
        args=(
            mgd.TempInputFile('infile_copy.obj'),
            mgd.TempInputFile('hmmcopy_res.obj'),
            mgd.TempSpace('correction_plots_dir'),
            mgd.TempSpace('hmmcopy_plots_dir'),
            mgd.OutputFile(bias_plots_pdf),
            mgd.OutputFile(correction_plots_pdf),
            mgd.OutputFile(hmmcopy_plots_pdf),
        )
    )

    workflow.transform(
        name='annot_hmm',
        func=tasks.annot_hmm,
        axes=('sample_id',),
        ctx={
            'mem': config['memory']['med'],
            'pool_id': config['pools']['standard'],
            'ncpus': 1,
        },
        args=(
            mgd.TempInputFile('hmmcopy_segments.txt'),
            mgd.OutputFile(pygene_outfile),
            config,
        )
    )

    return workflow
