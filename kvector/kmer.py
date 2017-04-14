from __future__ import print_function


import itertools

from Bio import SeqIO
import joblib
import numpy as np
import pandas as pd
import pybedtools

DIRECTIONS = 'upstream', 'downstream'


RNA = 'ACGU'
DNA = 'ACGT'


def interval_to_str_name(interval):
    """Convert a pybedtools interval to a string representing its position"""
    str_name = '{chrom}:{start}-{stop}'.format(
        chrom=interval.chrom, start=interval.start, stop=interval.stop)

    if interval.strand != '.':
        str_name += '({})'.format(interval.strand)
    return str_name


def score_kmers(pwm, kmers):
    """Generator to score kmers given a position-weight matrix

    Parameters
    ----------
    pwm : pandas.DataFrame
        A (length, 4) dataframe of the weight of each position's probability
        of each nucleotide
    kmers : list of str
        A list of kmers strings, e.g. ['AAA', 'AAC', 'AAG']

    """
    motif_length = pwm.shape[0]

    for kmer in map(list, kmers):
        k = len(kmer)

        divisor = min(k, motif_length)
        if k == motif_length:
            score = np.sum(pwm.lookup(range(motif_length), kmer))/divisor
        elif k > motif_length:
            starts = range(k - motif_length + 1)
            n_positions = len(starts)
            score = sum(np.sum(
                pwm.lookup(range(motif_length),
                           kmer[start:(start+motif_length)]))/divisor
                        for start in starts)/n_positions
        else:
            # k < motif_length
            starts = range(motif_length - k + 1)
            n_positions = len(starts)
            score = sum(np.sum(pwm.lookup(range(start, start+k), kmer))/divisor
                        for start in starts)/n_positions
        yield score


def count_kmers(filename, kmer_lengths=(4, 5, 6), format='fasta',
                residues=DNA):
    """Observe the number of substrings of specific lengths in sequence file

    A k-mer is a DNA (or RNA!) "word" of a specific length "k".

    Parameters
    ----------
    filename : str
        Name of a sequence file
    kmer_lengths : list of ints
        Lengths of the kmers to count. 8 or more (4^8 sequences) is generally
        not storable in memeory
    format : str, optional
        Format of the sequence file. Default is "fasta"
    residues : str, optional
        The residues to count, default is "ACGT" (DNA)

    Returns
    -------
    kmer_matrix : pandas.DataFrame
        A (kmers, sequences) dataframe of the kmers observed
    """
    if isinstance(kmer_lengths, int):
        kmer_lengths = [kmer_lengths]
    kmers = make_kmers(kmer_lengths, residues=residues)

    with open(filename) as f:
        records = list(SeqIO.parse(f, format))
    kmer_matrix = pd.DataFrame(0, columns=kmers, dtype=int,
                               index=range(len(records)))

    for col, record in enumerate(records):
        for k in kmer_lengths:
            for i in range(len(record) - k + 1):
                kmer = str(record[i:(i+k)].seq).upper()
                kmer_matrix.loc[col, kmer] += 1
    return kmer_matrix


def make_kmers(kmer_lengths, residues=DNA):
    """Create all possible substrings of provided lengths

    Parameters
    ----------
    kmer_lengths : list of ints
        Lengths of substrings to create
    residues : str
        The strings to sample from

    Returns
    -------
    kmers : list of str
        All possible kmers of the lengths provided, given the residues
    """
    try:
        return list(itertools.chain(
            *[map(
                lambda x: ''.join(x), itertools.product(
                    residues, repeat=k)) for k
              in kmer_lengths]))
    except TypeError:
        # Kmer length is only one number
        return list(map(lambda x: ''.join(x), itertools.product(
            residues, repeat=kmer_lengths)))


def _count_kmers_single_interval(interval, genome_fasta, intersect,
                                 kmer_lengths, residues):
    """Internal function for multithreaded counting of kmers"""
    # print('interval!', interval.chrom, interval.start, interval.stop)
    minibed = pybedtools.BedTool([interval])
    if intersect is not None:
        minibed = minibed.intersect(intersect)
    seqs = minibed.sequence(fi=genome_fasta, s=True)
    k = count_kmers(seqs.seqfn, kmer_lengths=kmer_lengths,
                    residues=residues).sum()
    k.name = interval.name
    return k


def per_interval_kmers(bed, fasta, intersect=None,
                       kmer_lengths=(4, 5, 6), residues=DNA, threads=-1):
    """Create a matrix of k-mer observations for each genomic region

    Parameters
    ----------
    bed : str or pybedtools.BedTool
        Either a filepath or pybedtools.BedTool of the genomic intervals whose
        kmers you want to count
    fasta : str
        Path to the genome fasta file
    intersect : str or pybedtools.BedTool
        Either a filepath or pybedtools.BedTool of another region location,
        e.g. conserved elements, that you want to intersect with when
        searching for k-mers
    threads : int
        Number of threads to use for multithreading of the job. Default is -1,
        which takes the maximum number of threads available. Usees joblib.

    Returns
    -------
    kmers : pandas.DataFrame
        A (n_intervals, n_kmers) matrix of the number of DNA words observed in
        each interval, possibly filtered on only the regions that intersect
        with the original intervals
    """
    if not isinstance(bed, pybedtools.BedTool):
        bed = pybedtools.BedTool(bed)

    if intersect is not None and not isinstance(intersect, pybedtools.BedTool):
        intersect = pybedtools.BedTool(intersect)

    if threads != 0:
        counts = joblib.Parallel(n_jobs=threads)(
            joblib.delayed(_count_kmers_single_interval)(
                interval, fasta, intersect, kmer_lengths, residues)
            for interval in bed)
    else:
        counts = []
        for interval in bed:
            counts.append(_count_kmers_single_interval(
                interval, fasta, intersect, kmer_lengths, residues))
    kmers = pd.concat(counts, axis=1).T
    index = [name if name != '.' else interval_to_str_name(interval)
             for name, interval in zip(kmers.index, bed)]
    kmers.index = index
    kmers = kmers.astype(int)
    return kmers
