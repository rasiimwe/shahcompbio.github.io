'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

import csv
import os
import random
import yaml

from pyclone.config import get_mutation
from pyclone.pyclone_beta_binomial import run_pyclone_beta_binomial_analysis
from pyclone.pyclone_binomial import run_pyclone_binomial_analysis
from pyclone.utils import make_directory, make_parent_directory

import pyclone.paths as paths
import pyclone.post_process as post_process
import pyclone.post_process.plot as plot

#=======================================================================================================================
# PyClone analysis
#=======================================================================================================================


def run_analysis_pipeline(args):
    config_file = _setup_analysis(
        density=args.density,
        in_files=args.in_files,
        init_method=args.init_method,
        num_iters=args.num_iters,
        samples=args.samples,
        prior=args.prior,
        tumour_contents=args.tumour_contents,
        working_dir=args.working_dir,
        config_extras_file=args.config_extras_file,
    )

    _run_analysis(config_file, args.seed)

    tables_dir = os.path.join(args.working_dir, 'tables')

    make_directory(tables_dir)

    for table_type in ['cluster', 'loci']:
        out_file = os.path.join(tables_dir, '{0}.tsv'.format(table_type))

        _build_table(
            config_file=config_file,
            out_file=out_file,
            burnin=args.burnin,
            max_clusters=args.max_clusters,
            mesh_size=args.mesh_size,
            table_type=table_type,
            thin=args.thin
        )

    plots_dir = os.path.join(args.working_dir, 'plots')

    plots = [
        ('cluster', 'density'),
        ('cluster', 'parallel_coordinates'),
        ('cluster', 'scatter'),
        ('loci', 'density'),
        ('loci', 'parallel_coordinates'),
        ('loci', 'scatter'),
        ('loci', 'similarity_matrix'),
        ('loci', 'vaf_parallel_coordinates'),
        ('loci', 'vaf_scatter')
    ]

    for category, plot_type in plots:

        plot_file = os.path.join(plots_dir, category, '{0}.{1}'.format(plot_type, args.plot_file_format))

        make_parent_directory(plot_file)

        if category == 'cluster':

            _cluster_plot(
                config_file,
                plot_file,
                args.burnin,
                args.max_clusters,
                args.mesh_size,
                args.min_cluster_size,
                plot_type,
                args.samples,
                args.thin
            )

        elif category == 'loci':

            _loci_plot(
                config_file,
                plot_file,
                plot_type,
                burnin=args.burnin,
                min_cluster_size=args.min_cluster_size,
                samples=args.samples,
                thin=args.thin
            )


def _write_config_file(
        config_file,
        density,
        init_method,
        mutations_files,
        num_iters,
        tumour_contents,
        working_dir,
        config_extras_file=None):

    config = {}

    config['num_iters'] = num_iters

    config['base_measure_params'] = {'alpha': 1, 'beta': 1}

    config['concentration'] = {
        'value': 1.0,
        'prior': {
            'shape': 1.0,
            'rate': 0.001
        }
    }

    config['density'] = density

    if density == 'pyclone_beta_binomial':
        config['beta_binomial_precision_params'] = {
            'value': 1000,
            'prior': {
                'shape': 1.0,
                'rate': 0.001
            },
            'proposal': {'precision': 0.01}
        }

    config['init_method'] = init_method

    config['working_dir'] = os.path.abspath(working_dir)

    config['trace_dir'] = 'trace'

    config['samples'] = {}

    for sample_id in mutations_files:
        config['samples'][sample_id] = {
            'mutations_file': mutations_files[sample_id],
            'tumour_content': {
                'value': tumour_contents[sample_id]
            },
            'error_rate': 0.001
        }

    if config_extras_file is not None:
        config.update(yaml.load(open(config_extras_file)))

    with open(config_file, 'w') as fh:
        yaml.dump(config, fh, default_flow_style=False, Dumper=Dumper)


def run_analysis(args):
    _run_analysis(args.config_file, args.seed)


def _run_analysis(config_file, seed):
    if seed is not None:
        random.seed(seed)

    config = paths.load_config(config_file)

    alpha = config['concentration']['value']

    if 'prior' in config['concentration']:
        alpha_priors = config['concentration']['prior']
    else:
        alpha_priors = None

    num_iters = config['num_iters']

    density = config['density']

    if density == 'pyclone_beta_binomial':
        run_pyclone_beta_binomial_analysis(
            config_file,
            num_iters,
            alpha,
            alpha_priors
        )

    elif density == 'pyclone_binomial':
        run_pyclone_binomial_analysis(
            config_file,
            num_iters,
            alpha,
            alpha_priors
        )

    else:
        raise Exception('{0} is not a valid density for PyClone.'.format(density))


def setup_analysis(args):
    _setup_analysis(
        config_extras_file=args.config_extras_file,
        density=args.density,
        in_files=args.in_files,
        init_method=args.init_method,
        num_iters=args.num_iters,
        samples=args.samples,
        prior=args.prior,
        tumour_contents=args.tumour_contents,
        working_dir=args.working_dir,
    )


def _setup_analysis(
        config_extras_file,
        density,
        in_files,
        init_method,
        num_iters,
        samples, prior,
        tumour_contents,
        working_dir):

    make_directory(working_dir)

    make_directory(os.path.join(working_dir, 'yaml'))

    mutations_files = OrderedDict()

    _tumour_contents = {}

    for i, in_file in enumerate(in_files):
        if samples is not None:
            sample_id = samples[i]

        else:
            sample_id = os.path.splitext(os.path.basename(in_file))[0]

        mutations_files[sample_id] = os.path.join(working_dir, 'yaml', '{0}.yaml'.format(sample_id))

        _build_mutations_file(
            in_file,
            mutations_files[sample_id],
            prior
        )

        if tumour_contents is not None:
            _tumour_contents[sample_id] = tumour_contents[i]

        else:
            _tumour_contents[sample_id] = 1.0

    config_file = os.path.join(working_dir, 'config.yaml')

    _write_config_file(
        config_file=config_file,
        density=density,
        init_method=init_method,
        mutations_files=mutations_files,
        num_iters=num_iters,
        tumour_contents=_tumour_contents,
        working_dir=working_dir,
        config_extras_file=config_extras_file,
    )

    return config_file

#=======================================================================================================================
# Input file code
#=======================================================================================================================


def build_mutations_file(args):
    _build_mutations_file(
        args.in_file,
        args.out_file,
        args.prior
    )


def _build_mutations_file(in_file, out_file, prior):
    config = {}

    reader = csv.DictReader(open(in_file), delimiter='\t')

    config['mutations'] = []

    for row in reader:
        mutation_id = row['mutation_id']

        ref_counts = int(row['ref_counts'])

        var_counts = int(row['var_counts'])

        normal_cn = int(row['normal_cn'])

        minor_cn = int(row['minor_cn'])

        major_cn = int(row['major_cn'])

        mutation = get_mutation(
            mutation_id,
            ref_counts,
            var_counts,
            normal_cn,
            minor_cn,
            major_cn,
            prior
        )

        config['mutations'].append(mutation.to_dict())

    make_parent_directory(out_file)

    fh = open(out_file, 'w')

    yaml.dump(config, fh, Dumper=Dumper)

    fh.close()

#=======================================================================================================================
# Post processing code
#=======================================================================================================================


def build_table(args):
    _build_table(
        config_file=args.config_file,
        out_file=args.out_file,
        burnin=args.burnin,
        max_clusters=args.max_clusters,
        mesh_size=args.mesh_size,
        table_type=args.table_type,
        thin=args.thin
    )


def _build_table(config_file, out_file, burnin, max_clusters, mesh_size, table_type, thin):
    if table_type == 'cluster':
        df = post_process.clusters.load_summary_table(
            config_file,
            burnin=burnin,
            max_clusters=max_clusters,
            mesh_size=mesh_size,
            thin=thin,
        )

    elif table_type == 'loci':
        df = post_process.loci.load_table(
            config_file,
            burnin,
            thin,
            max_clusters=max_clusters,
            old_style=False
        )

    elif table_type == 'old_style':
        df = post_process.loci.load_table(
            config_file,
            burnin,
            thin,
            max_clusters=max_clusters,
            old_style=True
        )

    df.to_csv(out_file, index=False, sep='\t')


def cluster_plot(args):
    _cluster_plot(
        config_file=args.config_file,
        plot_file=args.plot_file,
        burnin=args.burnin,
        max_clusters=args.max_clusters,
        mesh_size=args.mesh_size,
        min_cluster_size=args.min_cluster_size,
        plot_type=args.plot_type,
        samples=args.samples,
        thin=args.thin,
    )


def _cluster_plot(config_file, plot_file, burnin, max_clusters, mesh_size, min_cluster_size, plot_type, samples, thin):

    if plot_type == 'density':

        plot.clusters.density_plot(
            config_file,
            plot_file,
            burnin=burnin,
            thin=thin,
            max_clusters=max_clusters,
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples,
        )

    elif plot_type == 'parallel_coordinates':

        plot.clusters.parallel_coordinates_plot(
            config_file,
            plot_file,
            burnin=burnin,
            max_clusters=max_clusters,
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples,
            thin=thin
        )

    elif plot_type == 'scatter':

        plot.clusters.scatter_plot(
            config_file,
            plot_file,
            burnin=burnin,
            max_clusters=max_clusters,
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples,
            thin=thin
        )


def loci_plot(args):
    _loci_plot(
        args.config_file,
        args.plot_file,
        args.plot_type,
        burnin=args.burnin,
        max_clusters=args.max_clusters,
        min_cluster_size=args.min_cluster_size,
        samples=args.samples,
        thin=args.thin)


def _loci_plot(
        config_file,
        plot_file,
        plot_type,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1):

    kwargs = {
        'burnin': burnin,
        'max_clusters': max_clusters,
        'min_cluster_size': min_cluster_size,
        'samples': samples,
        'thin': thin
    }

    if plot_type.startswith('vaf'):

        kwargs['value'] = 'variant_allele_frequency'

    if plot_type == 'density':
        [kwargs.pop(x) for x in kwargs.keys() if 'cluster' in x]

        plot.loci.density_plot(
            config_file,
            plot_file,
            **kwargs
        )

    elif plot_type.endswith('parallel_coordinates'):

        plot.loci.parallel_coordinates_plot(
            config_file,
            plot_file,
            **kwargs
        )

    elif plot_type.endswith('scatter'):

        plot.loci.scatter_plot(
            config_file,
            plot_file,
            **kwargs
        )

    elif plot_type == 'similarity_matrix':

        plot.loci.similarity_matrix_plot(
            config_file,
            plot_file,
            **kwargs
        )
