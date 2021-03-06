import os
import sys

import pypeliner
import pypeliner.managed as mgd
import yaml
from wgs.utils import helpers
from wgs.workflows import realignment


def realign_bams(samples, inputs, outputs, out_dir, config, config_globals,single_node=False):
    outputs = dict([(sampid, outputs[sampid])
                    for sampid in samples])
    inputs = dict([(sampid, inputs[sampid])
                   for sampid in samples])

    os.path.join(out_dir, 'input.yaml')

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='realign_bam_file',
        func=realignment.realign_bam_files,
        args=(
            mgd.InputFile("input.bam", "sample_id", axes_origin=[], fnames=inputs),
            mgd.OutputFile("output.bam", "sample_id", axes_origin=[], fnames=outputs),
            out_dir,
            config,
            config_globals,
            samples
        ),
        kwargs={'single_node': single_node}
    )

    return workflow


def realign_bam_workflow(args):
    config = helpers.load_yaml(args['config_file'])
    config_globals = config['globals']
    config = config['alignment']

    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow(ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs']))

    outdir = args['out_dir']
    meta_yaml = os.path.join(outdir, 'metadata.yaml')
    input_yaml_blob = os.path.join(outdir, 'input.yaml')

    config = helpers.load_yaml(args['config_file'])
    config = config['alignment']

    yamldata = yaml.safe_load(open(args['input_yaml']))

    samples = yamldata.keys()

    input_bams = {sample: yamldata[sample]['input'] for sample in samples}
    output_bams = {sample: yamldata[sample]['output'] for sample in samples}

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name="realign",
        func=realign_bams,
        ctx=helpers.get_default_ctx(),
        args=(
            samples,
            mgd.InputFile("input.bam", 'sample_id', fnames=input_bams,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile("realigned.bam", 'sample_id', fnames=output_bams,
                           extensions=['.bai'], axes_origin=[]),
            args["out_dir"],
            config,
            config_globals
        ),
        kwargs={'single_node': args['single_node']}
    )

    outputted_filenames = helpers.expand_list(output_bams, samples, 'sample_id')

    workflow.transform(
        name='generate_meta_files_results',
        func='wgs.utils.helpers.generate_and_upload_metadata',
        args=(
            sys.argv[0:],
            args["out_dir"],
            outputted_filenames,
            mgd.OutputFile(meta_yaml)
        ),
        kwargs={
            'input_yaml_data': helpers.load_yaml(args['input_yaml']),
            'input_yaml': mgd.OutputFile(input_yaml_blob),
            'metadata': {'type': 'realignment'}
        }
    )

    pyp.run(workflow)
