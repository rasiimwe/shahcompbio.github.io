import csv
import sys
import logging
import os
import argparse

import pypeliner
import pypeliner.managed as mgd

import wrappers
import utils.download

destruct_directory = os.environ.get('DESTRUCT_PACKAGE_DIRECTORY', None)
if destruct_directory is None:
    raise Exception('please set the $DESTRUCT_PACKAGE_DIRECTORY environment variable to the root of the destruct package')


if __name__ == '__main__':

    import destruct_test
    import create_breakpoint_simulation

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('simconfig',
                           help='Simulation configuration filename')

    argparser.add_argument('outdir',
                           help='Output directory')

    argparser.add_argument('--config',
                           help='Configuration filename')

    argparser.add_argument('--chromosomes', nargs='*', type=str, default=['20'],
                           help='Reference chromosomes')

    argparser.add_argument('--include_nonchromosomal',  action='store_true',
                           help='Include non chromosomal reference sequences')

    args = vars(argparser.parse_args())

    config = {}

    if args['config'] is not None:
        execfile(args['config'], {}, config)

    config.update(args)

    pyp = pypeliner.app.Pypeline([destruct_test, create_breakpoint_simulation], config)

    try:
        os.makedirs(args['outdir'])
    except OSError:
        pass

    ctx = {'mem':4}

    pyp.sch.transform('read_params', (), ctx,
        destruct_test.read_simulation_params,
        mgd.TempOutputObj('simulation.params'),
        mgd.InputFile(args['simconfig']))

    pyp.sch.setobj(mgd.TempOutputObj('chromosomes'), args['chromosomes'])
    pyp.sch.setobj(mgd.TempOutputObj('include_nonchromosomal'), args['include_nonchromosomal'])

    pyp.sch.transform('create_genome', (), ctx,
        destruct_test.create_genome,
        None,
        mgd.TempInputObj('chromosomes'),
        mgd.TempInputObj('include_nonchromosomal'),
        mgd.OutputFile(os.path.join(args['outdir'], 'genome.fasta')))

    pyp.sch.transform('create_sim', (), ctx,
        create_breakpoint_simulation.create,
        None,
        mgd.TempInputObj('simulation.params'),
        mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.fasta')),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.tsv')),
        mgd.TempOutputFile('concordant.1.fastq'),
        mgd.TempOutputFile('concordant.2.fastq'),
        mgd.TempOutputFile('discordant.1.fastq'),
        mgd.TempOutputFile('discordant.2.fastq'))

    pyp.sch.commandline('cat1', (), ctx, 'cat', mgd.TempInputFile('concordant.1.fastq'), mgd.TempInputFile('discordant.1.fastq'), '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.1.fastq')))
    pyp.sch.commandline('cat2', (), ctx, 'cat', mgd.TempInputFile('concordant.2.fastq'), mgd.TempInputFile('discordant.2.fastq'), '>', mgd.OutputFile(os.path.join(args['outdir'], 'simulated.2.fastq')))

    bwaalign_script = os.path.join(destruct_directory, 'scripts', 'bwaalign.py')

    pyp.sch.commandline('bwa_align', (), ctx, 
        sys.executable,
        bwaalign_script,
        mgd.InputFile(os.path.join(args['outdir'], 'genome.fasta')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.1.fastq')),
        mgd.InputFile(os.path.join(args['outdir'], 'simulated.2.fastq')),
        mgd.TempOutputFile('simulated.unsorted.bam'),
        '--tmp', mgd.TempFile('bwa_tmp'))

    pyp.sch.transform('samtools_sort_index', (), ctx,
        destruct_test.samtools_sort_index,
        None,
        mgd.TempInputFile('simulated.unsorted.bam'),
        mgd.OutputFile(os.path.join(args['outdir'], 'simulated.bam')))

    pyp.run()


