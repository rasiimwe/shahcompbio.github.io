'''
Created on Nov 30, 2015

@author: Andrew Roth
'''
import matplotlib.gridspec as gs
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import seaborn as sb

import pyclone.post_process as post_process

import defaults
import _scatter
import utils


def density_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        mesh_size=101,
        min_cluster_size=0,
        samples=None,
        thin=1):

    df = post_process.clusters.load_table(
        config_file,
        burnin=burnin,
        thin=thin,
        max_clusters=max_clusters,
        mesh_size=mesh_size,
        min_size=min_cluster_size
    )

    sizes = df[['cluster_id', 'size']].drop_duplicates().set_index('cluster_id').to_dict()['size']

    if samples is None:
        samples = sorted(df['sample_id'].unique())

    else:
        df = df[df['sample_id'].isin(samples)]

    num_samples = len(samples)

    clusters = df['cluster_id'].unique()

    postions = range(1, len(clusters) + 1)

    utils.setup_plot()

    width = 8

    height = 2 * num_samples + 1

    fig = pp.figure(figsize=(width, height))

    grid = gs.GridSpec(nrows=num_samples, ncols=1)

    colors = utils.get_clusters_color_map(pd.Series(clusters))

    for ax_index, sample_id in enumerate(samples):
        plot_df = df[df['sample_id'] == sample_id]

        plot_df = plot_df.drop(['sample_id', 'size'], axis=1).set_index('cluster_id')

        ax = fig.add_subplot(grid[ax_index])

        utils.setup_axes(ax)

        ax.annotate(
            sample_id,
            xy=(1.01, 0.5),
            xycoords='axes fraction',
            fontsize=defaults.axis_label_font_size
        )

        for i, (cluster_id, log_pdf) in enumerate(plot_df.iterrows()):
            pos = postions[i]

            y = log_pdf.index.astype(float)

            x = np.exp(log_pdf)

            x = (x / x.max()) * 0.3

            ax.fill_betweenx(y, pos - x, pos + x, color=colors[cluster_id], where=(x > 1e-6))

        ax.set_xticks(postions)

        if ax_index == (num_samples - 1):
            x_tick_labels = ['{0} (n={1})'.format(x, sizes[x]) for x in clusters]

            ax.set_xticklabels(
                x_tick_labels,
                rotation=90
            )

            ax.set_xlabel(defaults.cluster_label, fontsize=defaults.axis_label_font_size)

        else:
            ax.set_xticklabels([])

        utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)

        ax.set_ylim(defaults.cellular_prevalence_limits)

    if num_samples == 1:
        ax.set_ylabel(
            defaults.cellular_prevalence_label,
            fontsize=defaults.axis_label_font_size
        )

    else:
        fig.text(
            -0.01,
            0.5,
            defaults.cellular_prevalence_label,
            fontsize=defaults.axis_label_font_size,
            ha='center',
            rotation=90,
            va='center'
        )

    grid.tight_layout(fig)

    utils.save_figure(fig, plot_file)


def parallel_coordinates_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        mesh_size=101,
        min_cluster_size=0,
        samples=None,
        thin=1):

    utils.setup_plot()

    plot_df = post_process.clusters.load_summary_table(
        config_file,
        burnin=burnin,
        max_clusters=max_clusters,
        mesh_size=mesh_size,
        min_size=min_cluster_size,
        thin=thin,
    )

    if samples is None:
        samples = sorted(plot_df['sample_id'].unique())

    else:
        plot_df = plot_df[plot_df['sample_id'].isin(samples)]

    clusters = sorted(plot_df['cluster_id'].unique())

    plot_df['sample_index'] = plot_df['sample_id'].apply(lambda x: samples.index(x))

    plot_df = plot_df.sort_values(by='sample_index')

    grid = sb.FacetGrid(
        plot_df,
        hue='cluster_id',
        hue_order=clusters,
        palette='husl'
    )

    grid.map(
        pp.errorbar,
        'sample_index',
        'mean',
        'std',
        marker=defaults.line_plot_marker,
        markersize=defaults.line_plot_marker_size
    )

    ax = grid.ax

    utils.setup_axes(ax)

    fig = grid.fig

    # Legend
    box = ax.get_position()

    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cluster')

    # Axis formatting
    ax.set_xticks(sorted(plot_df['sample_index'].unique()))

    ax.set_xticklabels(samples)

    ax.set_xlabel(defaults.sample_label, fontsize=defaults.axis_label_font_size)

    ax.set_ylabel(defaults.cellular_prevalence_label, fontsize=defaults.axis_label_font_size)

    utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)

    # Plot limits
    ax.set_xlim(
        plot_df['sample_index'].min() - 0.1,
        plot_df['sample_index'].max() + 0.1
    )

    ax.set_ylim(*defaults.cellular_prevalence_limits)

    # Resize and save figure
    fig.set_size_inches(*utils.get_parallel_coordinates_figure_size(samples))

    utils.save_figure(fig, plot_file)


def scatter_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        mesh_size=101,
        min_cluster_size=0,
        samples=None,
        thin=1):

    utils.setup_plot()

    df = post_process.clusters.load_summary_table(
        config_file,
        burnin=burnin,
        max_clusters=max_clusters,
        mesh_size=mesh_size,
        min_size=min_cluster_size,
        thin=thin,
    )

    mean_df = df.pivot(index='cluster_id', columns='sample_id', values='mean')

    error_df = df.pivot(index='cluster_id', columns='sample_id', values='std')

    if samples is None:
        samples = sorted(df['sample_id'].unique())

    color_map = utils.get_clusters_color_map(pd.Series(df['cluster_id']))

    _scatter.plot_all_pairs(
        color_map,
        mean_df,
        plot_file,
        samples,
        error_df=error_df
    )
