'''
Created on Feb 21, 2018

@author: pwalters
'''
import pypeliner
import pypeliner.managed as mgd
from wgs.utils import helpers


def create_annotation_workflow(
        input_vcf,
        annotated_vcf,
        global_config,
        varcall_config,
        vcftools_docker=None,
        snpeff_docker=None
):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='run_snpeff',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_snpeff',
        args=(
            mgd.InputFile(input_vcf),
            mgd.TempOutputFile('annotSnpEff.vcf'),
            varcall_config,
        ),
        kwargs={'docker_image': snpeff_docker}
    )

    workflow.transform(
        name='run_mutation_assessor',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['med'],
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_mutation_assessor',
        args=(
            mgd.TempInputFile('annotSnpEff.vcf'),
            mgd.TempOutputFile('annotMA.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_DBSNP',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_DBSNP',
        args=(
            mgd.TempInputFile('annotMA.vcf'),
            mgd.TempOutputFile('flagDBsnp.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_1000gen',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_1000gen',
        args=(
            mgd.TempInputFile('flagDBsnp.vcf'),
            mgd.TempOutputFile('flag1000gen.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='run_cosmic',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.workflows.vcf_annotation.tasks.run_cosmic',
        args=(
            mgd.TempInputFile('flag1000gen.vcf'),
            mgd.TempOutputFile('cosmic.vcf'),
            varcall_config,
        ),
    )

    workflow.transform(
        name='low_mappability_flag',
        func='wgs.workflows.vcf_annotation.tasks.flag_low_mappability',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        args=(
            mgd.TempInputFile('cosmic.vcf'),
            mgd.TempOutputFile('low_mapp.vcf'),
            varcall_config['mappability_ref']
        ),
    ),

    workflow.transform(
        name='finalize',
        ctx=helpers.get_default_ctx(
            memory=global_config['memory']['high'],
            walltime='8:00', ),
        func='wgs.utils.vcf_tasks.finalise_vcf',
        args=(
            mgd.TempInputFile('low_mapp.vcf'),
            mgd.OutputFile(annotated_vcf, extensions=['.csi', '.tbi']),
        ),
        kwargs={'docker_image': vcftools_docker}
    )

    return workflow
