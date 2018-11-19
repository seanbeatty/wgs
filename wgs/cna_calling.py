import os
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import titan


def remixt_workflow(tumour_path, normal_path, breakpoints, sample_id, remixt_refdata, outdir):
    workflow = pypeliner.workflow.Workflow()

    remixt_dir = os.path.join(outdir, 'remixt')
    remixt_config = {}

    remixt_results_filename = os.path.join(remixt_dir, 'results.h5')
    remixt_raw_dir = os.path.join(remixt_dir, 'raw_data')

    workflow.subworkflow(
        name='remixt',
        func="remixt.workflow.create_remixt_bam_workflow",
        args=(
            mgd.InputFile(breakpoints),
            {sample_id: mgd.InputFile(tumour_path),
             sample_id + 'N': mgd.InputFile(normal_path)},
            {sample_id: mgd.OutputFile(remixt_results_filename)},
            remixt_raw_dir,
            remixt_config,
            remixt_refdata,
        ),
        kwargs={
            'normal_id': sample_id + 'N',
        }
    )

    return workflow



def cna_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    samples = inputs.keys()
    tumours = {sample: inputs[sample]['tumour'] for sample in samples}
    normals = {sample: inputs[sample]['normal'] for sample in samples}
    breakpoints = {sample: inputs[sample]['breakpoints'] for sample in samples}

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='titan',
        func=titan.create_titan_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            args['out_dir'],
            config['globals'],
            config['cna_calling'],
            config['cna_calling']['titan_intervals'],
            mgd.InputInstance("sample_id")
        ),
    )

    workflow.subworkflow(
        name='remixt',
        func=remixt_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('normal_bam', 'sample_id', fnames=normals, extensions=['.bai']),
            mgd.InputFile('tumour_bam', 'sample_id', fnames=tumours, extensions=['.bai']),
            mgd.InputFile('breakpoints', 'sample_id', fnames=breakpoints),
            mgd.InputInstance('sample_id'),
            config['cna_calling']['remixt_refdata'],
            args['out_dir'],
        ),
    )



    pyp.run(workflow)