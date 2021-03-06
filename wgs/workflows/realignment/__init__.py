import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers
from wgs.workflows import alignment


def realign_bam_files(inputs, outputs, outdir, config, config_globals, samples, single_node=False):
    inputs = dict([(sample, inputs[sample]) for sample in samples])
    outputs = dict([(sample, outputs[sample]) for sample in samples])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples)

    workflow.transform(
        name='bam_to_fastq',
        ctx=helpers.get_default_ctx(
            walltime='72:00',
            disk=500
        ),
        func="wgs.workflows.realignment.tasks.split_by_rg",
        axes=('sample_id',),
        args=(
            mgd.InputFile('input.bam', 'sample_id', fnames=inputs),
            mgd.TempOutputFile("inputdata_read1.fastq.gz", 'sample_id', "readgroup"),
            mgd.TempOutputFile("inputdata_read2.fastq.gz", 'sample_id', "readgroup", axes_origin=[]),
            mgd.TempSpace("bamtofastq", 'sample_id')
        )
    )

    workflow.transform(
        name='get_sample_info',
        func="wgs.workflows.realignment.tasks.get_read_group",
        axes=('sample_id',),
        ret=mgd.TempOutputObj('sample_info', 'sample_id'),
        args=(
            mgd.InputFile('input.bam', 'sample_id', fnames=inputs),
        )
    )

    workflow.subworkflow(
        name='align_samples',
        func=alignment.align_samples,
        args=(
            config,
            config_globals,
            mgd.TempInputFile("inputdata_read1.fastq.gz", "sample_id", "readgroup", axes_origin=[]),
            mgd.TempInputFile("inputdata_read2.fastq.gz", "sample_id", "readgroup", axes_origin=[]),
            mgd.OutputFile('output.bam', 'sample_id', fnames=outputs, extensions=['.bai'], axes_origin=[]),
            outdir,
            mgd.TempInputObj('sample_info', 'sample_id', axes_origin=[])
        ),
        kwargs={'single_node': single_node}
    )

    return workflow
