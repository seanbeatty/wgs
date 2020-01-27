import os
import sys

import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import breakpoint_calling_consensus
from wgs.workflows import destruct_wgs
from wgs.workflows import lumpy


def sv_calling_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)

    config = helpers.load_yaml(args['config_file'])
    inputs = helpers.load_yaml(args['input_yaml'])

    meta_yaml = os.path.join(args["out_dir"], 'metadata.yaml')
    input_yaml_blob = os.path.join(args["out_dir"], 'input.yaml')

    tumours = helpers.get_values_from_input(inputs, 'tumour')
    normals = helpers.get_values_from_input(inputs, 'normal')
    samples = tumours.keys()

    sv_outdir = os.path.join(args['out_dir'], 'breakpoints', '{sample_id}')
    destruct_breakpoints = os.path.join(sv_outdir, 'destruct_breakpoints.csv.gz')
    destruct_library = os.path.join(sv_outdir, 'destruct_library.csv.gz')
    destruct_raw_breakpoints = os.path.join(sv_outdir, 'destruct_raw_breakpoints.csv.gz')
    destruct_raw_library = os.path.join(sv_outdir, 'destruct_raw_library.csv.gz')
    destruct_reads = os.path.join(sv_outdir, 'destruct_reads.csv.gz')
    lumpy_vcf = os.path.join(sv_outdir, 'lumpy.vcf')
    parsed_csv = os.path.join(sv_outdir, 'filtered_consensus_calls.csv.gz')

    single_node = args['single_node']

    global_config = config['globals']
    config = config['sv_calling']

    workflow = pypeliner.workflow.Workflow(
        ctx=helpers.get_default_ctx(docker_image=config['docker']['wgs'])
    )

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.subworkflow(
        name='destruct',
        func=destruct_wgs.create_destruct_wgs_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile("tumour.bam", 'sample_id', fnames=tumours,
                          extensions=['.bai'], axes_origin=[]),
            mgd.InputFile("normal.bam", 'sample_id', fnames=normals,
                          extensions=['.bai'], axes_origin=[]),
            mgd.OutputFile('destruct_raw_breakpoints', 'sample_id', template=destruct_raw_breakpoints),
            mgd.OutputFile('destruct_raw_library', 'sample_id', template=destruct_raw_library),
            mgd.OutputFile('destruct_breakpoints', 'sample_id', template=destruct_breakpoints),
            mgd.OutputFile('destruct_library', 'sample_id', template=destruct_library),
            mgd.OutputFile('destruct_reads', 'sample_id', template=destruct_reads),
            mgd.InputInstance('sample_id'),
            global_config,
            config
        ),
        kwargs={'single_node': single_node}
    )

    workflow.subworkflow(
        name='lumpy',
        func=lumpy.create_lumpy_workflow,
        axes=('sample_id',),
        args=(
            mgd.OutputFile('lumpy_vcf', 'sample_id', template=lumpy_vcf),
            global_config,
            config
        ),
        kwargs={
            'tumour_bam': mgd.InputFile(
                "tumour.bam", 'sample_id', fnames=tumours,
                extensions=['.bai'], axes_origin=[]),
            'normal_bam': mgd.InputFile(
                "normal.bam", 'sample_id', fnames=normals,
                extensions=['.bai'], axes_origin=[]),
            'single_node': single_node
        },
    )

    workflow.subworkflow(
        name="consensus_calling",
        func=breakpoint_calling_consensus.create_consensus_workflow,
        axes=('sample_id',),
        args=(
            mgd.InputFile('destruct_breakpoints', 'sample_id', template=destruct_breakpoints),
            mgd.InputFile('lumpy_vcf', 'sample_id', template=lumpy_vcf),
            mgd.OutputFile('consensus_calls', 'sample_id', template=parsed_csv),
            global_config,
            config,
            mgd.InputInstance('sample_id')
        ),
    )
    
    filenames = [destruct_breakpoints, destruct_library, destruct_raw_breakpoints,
        destruct_reads, lumpy_vcf, parsed_csv]

    outputted_filenames = helpers.expand_list(filenames, samples, "sample_id")

    workflow.transform(
        name='generate_meta_files_results',
        func=helpers.generate_and_upload_metadata,
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
