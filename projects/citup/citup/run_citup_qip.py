import argparse

import pypeliner
import pypeliner.managed as mgd

import citup.trees
import citup.tasks


if __name__ == '__main__':

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('input_freqs',
                           help='Input Mutation Frequencies')

    argparser.add_argument('input_clusters',
                           help='Input Mutation Clusters')

    argparser.add_argument('output_results',
                           help='Output Results HDF5 Store')

    argparser.add_argument('--min_nodes', type=int, default=1,
                           help='Output For All Trees')

    argparser.add_argument('--max_nodes', type=int, default=8,
                           help='Output For All Trees')

    argparser.add_argument('--max_children_per_node', type=int, default=100,
                           help='Output For All Trees')

    args = vars(argparser.parse_args())

    min_nodes = args['min_nodes']
    max_nodes = args['max_nodes']
    max_children_per_node = args['max_children_per_node']

    with open(args['input_clusters'], 'r') as f:
        num_clusters = len(set([a.rstrip() for a in f]))

    min_nodes = min(min_nodes, num_clusters + 1)
    max_nodes = min(max_nodes, num_clusters + 1)

    print 'min_nodes: {}, max_nodes: {}'.format(min_nodes, max_nodes)

    pyp = pypeliner.app.Pypeline(modules=[citup], config=args)

    workflow = pypeliner.workflow.Workflow(default_ctx={'mem': 4})

    workflow.transform(
        name='create_trees',
        func=citup.trees.create_trees,
        ret=mgd.TempOutputObj('trees', 'tree'),
        args=(
            min_nodes,
            max_nodes,
            max_children_per_node,
        ),
    )
    
    workflow.commandline(
        name='run_citup',
        axes=('tree',),
        ctx={'mem': 16},
        args=(
            'citupqip',
            mgd.TempInputObj('trees', 'tree').prop('unlabeled_tree_string'),
            mgd.InputFile(args['input_freqs']),
            mgd.InputFile(args['input_clusters']),
            mgd.TempOutputFile('results', 'tree'),
        ),
    )

    workflow.transform(
        name='select_optimal_tree',
        func=citup.tasks.select_optimal_tree,
        args=(
            mgd.InputFile(args['input_freqs']),
            mgd.TempInputObj('trees', 'tree'),
            mgd.TempInputFile('results', 'tree'),
            mgd.OutputFile(args['output_results']),
        ),
    )

    pyp.run(workflow)

