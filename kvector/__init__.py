# -*- coding: utf-8 -*-

__author__ = 'Olga Botvinnik'
__email__ = 'olga.botvinnik@gmail.com'
__version__ = '0.1.0'

from .io import read_motifs, create_metadata, motifs_to_kmer_vectors
from .kmer import count_kmers, score_kmers, make_kmers

__all__ = ['read_motifs', 'create_metadata', 'motifs_to_kmer_vectors',
           'count_kmers', 'score_kmers', 'make_kmers']
