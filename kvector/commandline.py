#!/usr/bin/env python
from __future__ import print_function

import click

import kvector

@click.command()
@click.argument('bed', type=click.Path(dir_okay=False))
@click.argument('fasta', type=click.Path(dir_okay=False))
@click.option('--intersect',
              help='Bed file of other regions, e.g. conserved elements, that '
                   'you want to intersect when searching for k-mers.')
@click.option('--kmer-lengths', default='4,5,6',
              help='How long of DNA words to search for (aka values of k). '
                   'Default is "4,5,6".')
@click.option('--residues', default=kvector.kmer.DNA,
              help="Which letters to search for in the fasta file. Default is "
                   "'{}'.".format(kvector.kmer.DNA))
@click.option('--threads', default=-1,
              help='Number of threads/processors to use for paralell processing'
                   ' of a multithreaded job. Default is -1, which uses the '
                   'maximum number of threads available, via the "joblib" '
                   'module.')
@click.version_option(version=kvector.__version__)
def cli(bed, fasta, intersect, kmer_lengths, residues, threads):
    """Counts k-mers in the bed intervals and writes a csv to stdout

    \b
    Parameters
    ----------
    bed : str
        Location of a bed file of the genomic intervals whose kmers you want
        to count
    fasta : str
        Path to the genome fasta file containing all chromosomes. Must be
        indexed (usually has a `.fai` file in the same directory, created
        using `faidx`).
    """

    kmer_lengths = map(int, kmer_lengths.split(','))

    kmers = kvector.per_interval_kmers(bed, fasta, intersect=intersect,
                                       kmer_lengths=kmer_lengths,
                                       residues=residues, threads=threads)
    click.echo(kmers.to_csv())


if __name__ == '__main__':
    cli()
