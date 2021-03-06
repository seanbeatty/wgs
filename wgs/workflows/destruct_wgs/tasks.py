import logging
import multiprocessing
import os

import pypeliner
import pypeliner.managed as mgd


def run_destruct_local(
        tempdir, tumour_bam, normal_bam,
        sample_id, raw_breakpoints, raw_library,
        reads, destruct_config, refdata_destruct,
        ncpus=None, docker_image=None
):
    pipelinedir = os.path.join(tempdir, 'pipeline')
    tmpdir = os.path.join(tempdir, 'tmp')

    if not ncpus:
        ncpus = multiprocessing.cpu_count()

    config = {'pipelinedir': pipelinedir, 'tmpdir': tmpdir,
              'submit': 'local', 'maxjobs': ncpus,
              'loglevel': 'DEBUG'}

    pyp = pypeliner.app.Pypeline(config=config)
    workflow = pypeliner.workflow.Workflow(ctx={'docker_image': docker_image})

    logging.getLogger().setLevel(logging.DEBUG)

    workflow.subworkflow(
        name='destruct_local_in_job',
        func='destruct.workflow.create_destruct_workflow',
        ctx={'docker_image': docker_image},
        args=(
            {sample_id: mgd.InputFile(tumour_bam),
             sample_id + 'N': mgd.InputFile(normal_bam)},
            mgd.OutputFile(raw_breakpoints),
            mgd.OutputFile(raw_library),
            mgd.OutputFile(reads),
            destruct_config,
            refdata_destruct
        )
    )

    pyp.run(workflow)
