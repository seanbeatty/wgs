'''
Created on Feb 19, 2018

@author: dgrewal
'''
import argparse
import json
import pypeliner
from wgs import __version__


def add_global_args(subparser):
    subparser.add_argument("--input_yaml",
                                 required=True,
                                 help='''yaml file with tumour, normal and sampleids''')

    subparser.add_argument("--out_dir",
                                required=True,
                                help='''Path to output directory.''')

    subparser.add_argument("--config_file",
                                help='''Path to the config file.''')

    subparser.add_argument("--config_override",
                                type=json.loads,
                                help='''json string to override the defaults in config''')

    pypeliner.app.add_arguments(subparser)

    return subparser





def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--version', action='version',
                        version='{version}'.format(version=__version__))

    subparsers = parser.add_subparsers()

    #================
    # variant calling
    #================
    variant_calling = subparsers.add_parser("variant_calling")
    variant_calling.set_defaults(which='variant_calling')
    variant_calling = add_global_args(variant_calling)


    #================
    # breakpoints calling
    #================
    sv_calling = subparsers.add_parser("sv_calling")
    sv_calling.set_defaults(which='sv_calling')
    sv_calling = add_global_args(sv_calling)

    #================
    # copy number calling
    #================
    cna_calling = subparsers.add_parser("cna_calling")
    cna_calling.set_defaults(which='cna_calling')
    cna_calling = add_global_args(cna_calling)

    #======================================
    # generates pipeline and batch configs
    #======================================
    generate_config = subparsers.add_parser("generate_config")
    generate_config.set_defaults(which='generate_config')

    generate_config.add_argument("--pipeline_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--batch_config",
                                 help='''output yaml file''')

    generate_config.add_argument("--config_override",
                                 type=json.loads,
                                 help='''json string to override the defaults in config''')


    args = vars(parser.parse_args())

    return args