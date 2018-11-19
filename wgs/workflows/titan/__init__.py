import os
import pypeliner
import pypeliner.managed as mgd
import tasks
import yaml

def create_titan_workflow(normal_bam, tumour_bam, outdir,
                          global_config, config, intervals, sample_id):

    igv_template = os.path.join(outdir, sample_id, 'clusters_{numclusters}', 'ploidy_{ploidy}', 'igv_segs.txt')
    outfile_template = os.path.join(outdir, sample_id, 'clusters_{numclusters}', 'ploidy_{ploidy}', 'titan_markers.txt')
    params_template = os.path.join(outdir, sample_id, 'clusters_{numclusters}', 'ploidy_{ploidy}', 'titan_params.txt')
    segs_template = os.path.join(outdir, sample_id, 'clusters_{numclusters}', 'ploidy_{ploidy}', 'titan_segs.txt')
    plots_template = os.path.join(outdir, sample_id, 'clusters_{numclusters}', 'ploidy_{ploidy}', 'titan_plots.tar.gz')

    chunks = [(v['num_clusters'], v['ploidy']) for v in intervals]

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('numclusters', 'ploidy'),
        value=chunks,
    )


    workflow.transform(
        name='generate_intervals',
        func=tasks.generate_intervals,
        ctx={'mem': global_config['memory']['low'],
             'pool_id': global_config['pools']['standard'],
             'ncpus': 1, 'walltime': '02:00'},
             # 'walltime_retry_increment': 2},
        ret=mgd.OutputChunks('interval'),
        args=(
            config['reference_genome'],
            config['chromosomes']
        )
    )

    workflow.transform(
        name='run_museq',
        ctx={'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['multicore'],
             'ncpus': global_config['threads'], 'walltime': '02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.run_museq,
        axes=('interval',),
        args=(
            mgd.InputFile(tumour_bam, extensions=['.bai']),
            mgd.InputFile(normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('museq.vcf', 'interval'),
            mgd.TempOutputFile('museq.log', 'interval'),
            config,
            mgd.InputInstance('interval')
        ),
    )


    workflow.transform(
        name='merge_vcfs',
        ctx={'num_retry': 3, 'mem_retry_increment': 2,
             'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['multicore'],
             'ncpus': 1},
             # 'walltime': '02:00', 'walltime_retry_increment': 2},
        func=tasks.merge_vcfs,
        args=(
            mgd.TempInputFile('museq.vcf', 'interval'),
            mgd.TempOutputFile('museq.vcf'),
        )
    )


    workflow.transform(
        name='convert_museq_vcf2counts',
        ctx={'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['highmem'],
             'ncpus': 1, 'walltime': '02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.convert_museq_vcf2counts,
        args=(
            mgd.TempInputFile('museq.vcf'),
            mgd.TempOutputFile('museq_postprocess.txt'),
            config,
        ),
    )

    workflow.transform(
        name='run_readcounter_tumour',
        ctx={'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['highmem'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.run_readcounter,
        args=(
            mgd.InputFile(tumour_bam, extensions=['.bai']),
            mgd.TempOutputFile('tumour.wig'),
            config,
        ),
    )

    workflow.transform(
        name='run_readcounter_normal',
        ctx={'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['highmem'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.run_readcounter,
        args=(
            mgd.InputFile(normal_bam, extensions=['.bai']),
            mgd.TempOutputFile('normal.wig'),
            config,
        ),
    )

    workflow.transform(
        name='calc_correctreads_wig',
        ctx={'mem': global_config['memory']['low'],
             'pool_id': global_config['pools']['standard'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.calc_correctreads_wig,
        args=(
            mgd.TempInputFile('tumour.wig'),
            mgd.TempInputFile('normal.wig'),
            mgd.TempOutputFile('correct_reads.txt'),
            config,
        ),
    )

    workflow.transform(
        name='run_titan',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['high'],
             'pool_id': global_config['pools']['highmem'],
             'ncpus': 1, 'walltime': '06:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.run_titan,
        args=(
            mgd.TempInputFile('museq_postprocess.txt'),
            mgd.TempInputFile('correct_reads.txt'),
            mgd.OutputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
            mgd.TempOutputFile('titan.Rdata','numclusters', 'ploidy'),
            mgd.OutputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
            config['titan_params'],
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy')
        )
    )

    workflow.transform(
        name='plot_titan',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['low'],
             'pool_id': global_config['pools']['standard'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.plot_titan,
        args=(
            mgd.TempInputFile('titan.Rdata','numclusters', 'ploidy'),
            mgd.InputFile('titan_params', 'numclusters', 'ploidy', template=params_template),
            mgd.OutputFile('titan_plots', 'numclusters', 'ploidy', template=plots_template),
            mgd.TempSpace("titan_plots_tempdir", 'numclusters', 'ploidy'),
            config,
            mgd.InputInstance('numclusters'),
            mgd.InputInstance('ploidy')
        ),
    )

    workflow.transform(
        name='calc_cnsegments_titan',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['low'],
             'pool_id': global_config['pools']['standard'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.calc_cnsegments_titan,
        args=(
            mgd.InputFile('titan_outfile', 'numclusters', 'ploidy', template=outfile_template),
            mgd.OutputFile('titan_igv', 'numclusters', 'ploidy', template=igv_template),
            mgd.TempOutputFile('segs.csv', 'numclusters', 'ploidy'),
        ),
    )
    
    workflow.transform(
        name='annot_pygenes',
        axes=('numclusters', 'ploidy'),
        ctx={'mem': global_config['memory']['low'],
             'pool_id': global_config['pools']['standard'],
             'ncpus': 1, 'walltime':'02:00'},
             # 'walltime_retry_increment': 2},
        func=tasks.annot_pygenes,
        args=(
            mgd.TempInputFile('segs.csv', 'numclusters', 'ploidy'),
            mgd.OutputFile('titan_segs.csv', 'numclusters', 'ploidy', template=segs_template),
            config,
        ),
    )

    return workflow