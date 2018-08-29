#!/usr/bin/env python

#=======================================================================================================================
# PyClone
# Author : Andrew Roth
#=======================================================================================================================
import pyclone.run as run

import argparse


def main():
    parser = argparse.ArgumentParser(prog='PyClone')

    parser.add_argument('--version', action='version', version='PyClone-0.13.1')

    subparsers = parser.add_subparsers()

#----------------------------------------------------------------------------------------------------------------------
    setup_analysis_parser = subparsers.add_parser(
        'setup_analysis',
        help='''Setup a config file and mutations files for a PyClone analysis.'''
    )

    _setup_setup_analysis_parser(setup_analysis_parser)

#----------------------------------------------------------------------------------------------------------------------
    run_analysis_parser = subparsers.add_parser(
        'run_analysis',
        help='''Run an MCMC sampler to sample from the posterior of the PyClone model.'''
    )

    _setup_run_analysis_parser(run_analysis_parser)

#----------------------------------------------------------------------------------------------------------------------
    analysis_pipeline_parser = subparsers.add_parser(
        'run_analysis_pipeline',
        help='''Run a full PyClone analysis.'''
    )

    _setup_analysis_pipeline_parser(analysis_pipeline_parser)

#----------------------------------------------------------------------------------------------------------------------
    build_prior_parser = subparsers.add_parser(
        'build_mutations_file',
        help='''Build a YAML format file with mutation data and states prior to be used for PyClone analysis.'''
    )

    _setup_build_prior_parser(build_prior_parser)


#----------------------------------------------------------------------------------------------------------------------
    plot_clusters_parser = subparsers.add_parser(
        'plot_clusters',
        help='''Plot features of the clusters.'''
    )

    _setup_cluster_plot_parser(plot_clusters_parser)

#----------------------------------------------------------------------------------------------------------------------
    plot_loci_parser = subparsers.add_parser(
        'plot_loci',
        help='''Plot features of the loci.'''
    )

    _setup_loci_plot_parser(plot_loci_parser)

#----------------------------------------------------------------------------------------------------------------------
    build_table_parser = subparsers.add_parser(
        'build_table',
        help='''Build results table which contains cluster ids and (mean) cellular prevalence estimates.''')

    _setup_build_table_parser(build_table_parser)

#----------------------------------------------------------------------------------------------------------------------
    args = parser.parse_args()

    args.func(args)


def _setup_run_analysis_parser(parser):

    _add_config_file_args(parser)

    _add_seed_args(parser)

    parser.set_defaults(func=run.run_analysis)


def _setup_setup_analysis_parser(parser):

    parser.add_argument(
        '--in_files',
        nargs='+',
        required=True,
        help='''Space delimited list of tsv format files with copy number and allele count information. See 
        build_mutations_file command for information.'''
    )

    parser.add_argument(
        '--working_dir',
        required=True,
        help='''Path of directory where analysis pipeline files will be placed.'''
    )

    parser.add_argument(
        '--tumour_contents',
        nargs='+',
        type=float,
        default=None,
        help='''Space delimited list of tumour contents. Should match the order of --in_files. If not given tumour
        content is assumed to 1.0 in all samples.'''
    )

    parser.add_argument(
        '--samples',
        nargs='+', default=None, help='''Space delimited list of sample names. Should be in the same order as
        --in_files. If not set sample name will be inferred from file names and ordering in plots will be arbitrary.''')

    parser.add_argument(
        '--config_extras_file',
        required=False,
        help='''Path to configuration file with extra parameters used for analysis.'''
    )

    parser.add_argument(
        '--density',
        choices=['pyclone_binomial', 'pyclone_beta_binomial'],
        default='pyclone_beta_binomial',
        help='''Emission density for the model. Default is pyclone_beta_binomial.'''
    )

    _add_init_method_args(parser)

    parser.add_argument(
        '--num_iters',
        default=10000,
        type=int,
        help='''Number of iterations of the MCMC sampler to perform. Default is 10,000.'''
    )

    _add_prior_args(parser)

    parser.set_defaults(func=run.setup_analysis)


def _setup_analysis_pipeline_parser(parser):

    _setup_setup_analysis_parser(parser)

    _add_post_process_args(parser)

    _add_seed_args(parser)

    parser.add_argument(
        '--plot_file_format',
        default='pdf',
        choices=['pdf', 'svg'],
        help='''File format for plots. Default is pdf.'''
    )

    _add_max_clusters_args(parser)

    _add_mesh_size_args(parser)

    _add_min_cluster_size_args(parser)

    parser.set_defaults(func=run.run_analysis_pipeline)


def _setup_build_prior_parser(parser):

    parser.add_argument(
        '--in_file',
        required=True,
        help='''Path to tab separated input file. The input file should have header and the following columns: 
        mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn. Any additional columns will be ignored. 
        See examples for format.'''
    )

    parser.add_argument(
        '--out_file',
        required=True,
        help='''Path where YAML formatted PyClone input file will be written.'''
    )

    _add_prior_args(parser)

    parser.set_defaults(func=run.build_mutations_file)


def _setup_build_table_parser(parser):

    _add_config_file_args(parser)

    parser.add_argument(
        '--out_file',
        required=True,
        help='''Path where table will be written in tsv format.'''
    )

    parser.add_argument(
        '--table_type',
        choices=['cluster', 'loci', 'old_style'],
        required=True,
        help='''Build a table of results. Choices are: `cluster` for cluster specific information; `loci` for loci 
        specific information; `old_style` matches the 0.12.x PyClone output.'''
    )

    _add_max_clusters_args(parser)

    _add_mesh_size_args(parser)

    _add_post_process_args(parser)

    parser.set_defaults(func=run.build_table)


def _setup_cluster_plot_parser(parser):

    _add_config_file_args(parser)

    _add_plot_out_file_args(parser)

    parser.add_argument(
        '--plot_type',
        choices=['density', 'parallel_coordinates', 'scatter'],
        default='density',
        required=True,
        help='''Determines which style of plot will be done.'''
    )

    _add_post_process_args(parser)

    _add_max_clusters_args(parser)

    _add_mesh_size_args(parser)

    _add_min_cluster_size_args(parser)

    parser.add_argument(
        '--samples',
        default=None,
        nargs='+',
        help='''Samples to plot and order to plot them in.'''
    )

    parser.set_defaults(func=run.cluster_plot)


def _setup_loci_plot_parser(parser):

    _add_config_file_args(parser)

    _add_plot_out_file_args(parser)

    parser.add_argument(
        '--plot_type',
        choices=[
            'density',
            'parallel_coordinates',
            'scatter',
            'similarity_matrix',
            'vaf_parallel_coordinates',
            'vaf_scatter'
        ],
        default='density',
        required=True,
        help='''Determines which style of plot will be done.'''
    )

    _add_post_process_args(parser)

    _add_min_cluster_size_args(parser)

    _add_max_clusters_args(parser)

    parser.add_argument(
        '--samples',
        default=None,
        nargs='+',
        help='''Samples to plot and order to plot them in.'''
    )

    parser.set_defaults(func=run.loci_plot)

#=======================================================================================================================
# Common args
#=======================================================================================================================


def _add_config_file_args(parser):

    parser.add_argument(
        '--config_file',
        required=True,
        help='''Path to configuration file used for analysis.'''
    )


def _add_plot_out_file_args(parser):

    parser.add_argument(
        '--plot_file',
        required=True,
        help='Path to file where plot will be saved. Format can be controlled by changing file extension.'
    )


def _add_seed_args(parser):

    parser.add_argument(
        '--seed',
        default=None,
        type=int,
        help='''Set random seed so results can be reproduced. By default a random seed is chosen.'''
    )


def _add_prior_args(parser):

    parser.add_argument(
        '--prior',
        choices=['major_copy_number', 'parental_copy_number', 'total_copy_number'],
        default='major_copy_number',
        help='''Method used to set the possible genotypes. See online help for description. Default is major_copy_number.'''
    )


def _add_post_process_args(parser):

    parser.add_argument(
        '--burnin',
        default=0,
        type=int,
        help='''Number of samples to discard as burning for the MCMC chain. Default is 0.'''
    )

    parser.add_argument(
        '--thin',
        default=1,
        type=int,
        help='''Number of samples to thin MCMC trace. For example if thin=10 every tenth sample after burning will be 
        used for inference. Default is 1.'''
    )


def _add_min_cluster_size_args(parser):
    parser.add_argument(
        '--min_cluster_size',
        default=0,
        type=int,
        help='''Clusters with fewer mutations than this value will not be plotted.'''
    )


def _add_mesh_size_args(parser):
    parser.add_argument(
        '--mesh_size',
        default=101,
        type=int,
        help='''Number of points to use for approximating the cluster posteriors. Default is 101.'''
    )


def _add_init_method_args(parser):
    parser.add_argument(
        '--init_method', choices=['connected', 'disconnected'], default='disconnected',
        help='''How to initialise the DP clustering algorithm. `connected` places all data points in one cluster,
        preferred for large datasets. `disconnected places each data point in a separate cluster. Default
        `disconnected`.'''
    )


def _add_max_clusters_args(parser):
    parser.add_argument(
        '--max_clusters', default=None, type=int,
        help='''Maximum number of clusters to consider for post-processing. Note this does not affect the DP sampling
        only the final post-processing steps to get hard cluster assignments.'''
    )


if __name__ == '__main__':
    main()
