import itertools

import numpy as np

RNA = 'ACGU'
DNA = 'ACGT'


def score_kmers(pwm, kmers):
    """Generator to score kmers given a position-weight matrix

    Parameters
    ----------
    pwm : pandas.DataFrame
        A (length, 4) dataframe of the weight of each position's probability
        of each nucleotide
    kmers : list of list
        A list of kmers strings as lists, e.g. [['G', 'G', 'G', 'G', 'G', 'G']]

    """
    motif_length = pwm.shape[0]
    for kmer in kmers:
        k = len(kmer)

        divisor = min(k, motif_length)
        if k == motif_length:
            score = np.sum(pwm.lookup(range(motif_length), kmer))/divisor
        elif k > motif_length:
            starts = range(k - motif_length + 1)
            n_positions = len(starts)
            score = sum(np.sum(pwm.lookup(range(motif_length),
                                                kmer[start:(start+motif_length)]))/divisor
                        for start in starts)/n_positions
        else:
            # k < motif_length
            starts = range(motif_length - k + 1)
            n_positions = len(starts)
            score = sum(np.sum(pwm.lookup(range(start, start+k), kmer))/divisor
                        for start in starts)/n_positions
        yield score

def make_kmers(kmer_lengths, residues=DNA):
    return list(itertools.chain(
        *[map(lambda x: ''.join(x), itertools.product(residues, repeat=k)) for k
          in kmer_lengths]))