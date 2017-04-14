#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_kmer
----------------------------------

Tests for `kvector.kmer` module.
"""
import os

import pandas as pd
import pandas.util.testing as pdt
import pybedtools
import pytest
import six


@pytest.fixture
def fasta(data_folder):
    return '{}/example.fasta'.format(data_folder)


@pytest.fixture
def pwm():
    s = ''',A,C,G,T
    0,0.39532879396435,0.105513888686126,0.105513888686126,0.39364342774540506
    1,0.00770456803068082,0.00770456803068082,0.00770456803068082,0.976886297348457
    2,0.976886297348457,0.00770456803068082,0.00770456803068082,0.00770456803068082
    3,0.976886297348457,0.00770456803068082,0.00770456803068082,0.00770456803068082
    4,0.00770456803068082,0.00770456803068082,0.00770456803068082,0.976886297348457
    5,0.00770456803068082,0.00770456803068082,0.00770456803068082,0.976886297348457
    6,0.321131137484576,0.14380369811499302,0.478370946765367,0.0566942181612342
    '''
    return pd.read_csv(six.StringIO(s), index_col=0)

def test_interval_to_str_name_stranded():
    import kvector

    chrom, start, stop, strand = 'beyonce', 100, 200, '+'
    interval = pybedtools.Interval(chrom, start, stop,
                                   strand=strand)
    test = kvector.kmer.interval_to_str_name(interval)
    true = 'beyonce:100-200(+)'

    assert test == true


def test_interval_to_str_name_unstranded():
    import kvector

    chrom, start, stop= 'beyonce', 100, 200
    interval = pybedtools.Interval(chrom, start, stop)
    test = kvector.kmer.interval_to_str_name(interval)
    true = 'beyonce:100-200'

    assert test == true


def test_make_kmers(kmer_lengths):
    from kvector import make_kmers
    test = make_kmers(kmer_lengths)

    if kmer_lengths == 3:
        true = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA',
                'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC',
                'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG',
                'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT',
                'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
                'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC',
                'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG',
                'TTT']
    else:
        true = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA',
                'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC',
                'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG',
                'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT',
                'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
                'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC',
                'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG',
                'TTT', 'AAAA', 'AAAC', 'AAAG', 'AAAT', 'AACA', 'AACC', 'AACG',
                'AACT', 'AAGA', 'AAGC', 'AAGG', 'AAGT', 'AATA', 'AATC', 'AATG',
                'AATT', 'ACAA', 'ACAC', 'ACAG', 'ACAT', 'ACCA', 'ACCC', 'ACCG',
                'ACCT', 'ACGA', 'ACGC', 'ACGG', 'ACGT', 'ACTA', 'ACTC', 'ACTG',
                'ACTT', 'AGAA', 'AGAC', 'AGAG', 'AGAT', 'AGCA', 'AGCC', 'AGCG',
                'AGCT', 'AGGA', 'AGGC', 'AGGG', 'AGGT', 'AGTA', 'AGTC', 'AGTG',
                'AGTT', 'ATAA', 'ATAC', 'ATAG', 'ATAT', 'ATCA', 'ATCC', 'ATCG',
                'ATCT', 'ATGA', 'ATGC', 'ATGG', 'ATGT', 'ATTA', 'ATTC', 'ATTG',
                'ATTT', 'CAAA', 'CAAC', 'CAAG', 'CAAT', 'CACA', 'CACC', 'CACG',
                'CACT', 'CAGA', 'CAGC', 'CAGG', 'CAGT', 'CATA', 'CATC', 'CATG',
                'CATT', 'CCAA', 'CCAC', 'CCAG', 'CCAT', 'CCCA', 'CCCC', 'CCCG',
                'CCCT', 'CCGA', 'CCGC', 'CCGG', 'CCGT', 'CCTA', 'CCTC', 'CCTG',
                'CCTT', 'CGAA', 'CGAC', 'CGAG', 'CGAT', 'CGCA', 'CGCC', 'CGCG',
                'CGCT', 'CGGA', 'CGGC', 'CGGG', 'CGGT', 'CGTA', 'CGTC', 'CGTG',
                'CGTT', 'CTAA', 'CTAC', 'CTAG', 'CTAT', 'CTCA', 'CTCC', 'CTCG',
                'CTCT', 'CTGA', 'CTGC', 'CTGG', 'CTGT', 'CTTA', 'CTTC', 'CTTG',
                'CTTT', 'GAAA', 'GAAC', 'GAAG', 'GAAT', 'GACA', 'GACC', 'GACG',
                'GACT', 'GAGA', 'GAGC', 'GAGG', 'GAGT', 'GATA', 'GATC', 'GATG',
                'GATT', 'GCAA', 'GCAC', 'GCAG', 'GCAT', 'GCCA', 'GCCC', 'GCCG',
                'GCCT', 'GCGA', 'GCGC', 'GCGG', 'GCGT', 'GCTA', 'GCTC', 'GCTG',
                'GCTT', 'GGAA', 'GGAC', 'GGAG', 'GGAT', 'GGCA', 'GGCC', 'GGCG',
                'GGCT', 'GGGA', 'GGGC', 'GGGG', 'GGGT', 'GGTA', 'GGTC', 'GGTG',
                'GGTT', 'GTAA', 'GTAC', 'GTAG', 'GTAT', 'GTCA', 'GTCC', 'GTCG',
                'GTCT', 'GTGA', 'GTGC', 'GTGG', 'GTGT', 'GTTA', 'GTTC', 'GTTG',
                'GTTT', 'TAAA', 'TAAC', 'TAAG', 'TAAT', 'TACA', 'TACC', 'TACG',
                'TACT', 'TAGA', 'TAGC', 'TAGG', 'TAGT', 'TATA', 'TATC', 'TATG',
                'TATT', 'TCAA', 'TCAC', 'TCAG', 'TCAT', 'TCCA', 'TCCC', 'TCCG',
                'TCCT', 'TCGA', 'TCGC', 'TCGG', 'TCGT', 'TCTA', 'TCTC', 'TCTG',
                'TCTT', 'TGAA', 'TGAC', 'TGAG', 'TGAT', 'TGCA', 'TGCC', 'TGCG',
                'TGCT', 'TGGA', 'TGGC', 'TGGG', 'TGGT', 'TGTA', 'TGTC', 'TGTG',
                'TGTT', 'TTAA', 'TTAC', 'TTAG', 'TTAT', 'TTCA', 'TTCC', 'TTCG',
                'TTCT', 'TTGA', 'TTGC', 'TTGG', 'TTGT', 'TTTA', 'TTTC', 'TTTG',
                'TTTT']

    pdt.assert_equal(test, true)


def test_count_kmers(fasta, kmer_lengths):
    from kvector.kmer import count_kmers

    test = count_kmers(fasta, kmer_lengths=kmer_lengths)

    if kmer_lengths == 3:
        s = ''',AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT
0,2,3,1,0,0,2,1,2,1,0,3,0,0,0,0,0,1,1,1,0,1,10,2,5,0,0,3,1,2,5,1,0,1,1,2,0,1,3,0,0,3,3,5,0,0,2,0,0,2,0,0,0,1,4,1,1,0,0,0,1,0,0,0,0
1,2,3,1,0,0,2,1,2,1,0,3,0,0,0,0,0,1,1,1,0,1,12,3,5,0,1,3,1,2,5,1,0,1,1,2,0,1,4,0,0,3,3,5,0,0,2,0,0,2,0,0,0,1,4,1,1,0,0,0,1,0,0,0,0
2,6,5,6,0,0,3,2,4,3,5,3,1,0,0,0,0,2,2,2,0,3,16,7,7,3,3,6,4,2,8,4,0,7,2,4,0,2,8,3,2,7,5,7,0,0,3,1,2,2,0,0,0,2,6,4,1,0,1,3,1,0,2,0,0
3,6,6,4,1,12,2,0,5,3,9,2,3,1,2,5,0,5,12,8,4,4,12,2,6,0,2,2,0,2,8,8,2,4,0,4,3,10,7,2,5,4,3,3,3,1,0,4,2,2,1,1,0,3,3,0,4,4,10,6,1,0,1,3,2
4,19,9,1,7,7,8,0,8,4,1,1,3,2,6,5,6,9,8,4,6,4,4,0,8,0,1,2,0,3,8,6,10,3,2,2,1,7,0,3,3,2,2,0,5,4,1,7,3,4,4,2,5,9,5,0,8,2,9,6,7,6,7,6,18
'''
    else:
        s = ''',AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT,AAAA,AAAC,AAAG,AAAT,AACA,AACC,AACG,AACT,AAGA,AAGC,AAGG,AAGT,AATA,AATC,AATG,AATT,ACAA,ACAC,ACAG,ACAT,ACCA,ACCC,ACCG,ACCT,ACGA,ACGC,ACGG,ACGT,ACTA,ACTC,ACTG,ACTT,AGAA,AGAC,AGAG,AGAT,AGCA,AGCC,AGCG,AGCT,AGGA,AGGC,AGGG,AGGT,AGTA,AGTC,AGTG,AGTT,ATAA,ATAC,ATAG,ATAT,ATCA,ATCC,ATCG,ATCT,ATGA,ATGC,ATGG,ATGT,ATTA,ATTC,ATTG,ATTT,CAAA,CAAC,CAAG,CAAT,CACA,CACC,CACG,CACT,CAGA,CAGC,CAGG,CAGT,CATA,CATC,CATG,CATT,CCAA,CCAC,CCAG,CCAT,CCCA,CCCC,CCCG,CCCT,CCGA,CCGC,CCGG,CCGT,CCTA,CCTC,CCTG,CCTT,CGAA,CGAC,CGAG,CGAT,CGCA,CGCC,CGCG,CGCT,CGGA,CGGC,CGGG,CGGT,CGTA,CGTC,CGTG,CGTT,CTAA,CTAC,CTAG,CTAT,CTCA,CTCC,CTCG,CTCT,CTGA,CTGC,CTGG,CTGT,CTTA,CTTC,CTTG,CTTT,GAAA,GAAC,GAAG,GAAT,GACA,GACC,GACG,GACT,GAGA,GAGC,GAGG,GAGT,GATA,GATC,GATG,GATT,GCAA,GCAC,GCAG,GCAT,GCCA,GCCC,GCCG,GCCT,GCGA,GCGC,GCGG,GCGT,GCTA,GCTC,GCTG,GCTT,GGAA,GGAC,GGAG,GGAT,GGCA,GGCC,GGCG,GGCT,GGGA,GGGC,GGGG,GGGT,GGTA,GGTC,GGTG,GGTT,GTAA,GTAC,GTAG,GTAT,GTCA,GTCC,GTCG,GTCT,GTGA,GTGC,GTGG,GTGT,GTTA,GTTC,GTTG,GTTT,TAAA,TAAC,TAAG,TAAT,TACA,TACC,TACG,TACT,TAGA,TAGC,TAGG,TAGT,TATA,TATC,TATG,TATT,TCAA,TCAC,TCAG,TCAT,TCCA,TCCC,TCCG,TCCT,TCGA,TCGC,TCGG,TCGT,TCTA,TCTC,TCTG,TCTT,TGAA,TGAC,TGAG,TGAT,TGCA,TGCC,TGCG,TGCT,TGGA,TGGC,TGGG,TGGT,TGTA,TGTC,TGTG,TGTT,TTAA,TTAC,TTAG,TTAT,TTCA,TTCC,TTCG,TTCT,TTGA,TTGC,TTGG,TTGT,TTTA,TTTC,TTTG,TTTT
0,2,3,1,0,0,2,1,2,1,0,3,0,0,0,0,0,1,1,1,0,1,10,2,5,0,0,3,1,2,5,1,0,1,1,2,0,1,3,0,0,3,3,5,0,0,2,0,0,2,0,0,0,1,4,1,1,0,0,0,1,0,0,0,0,0,2,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,3,2,4,0,0,1,1,1,3,1,0,0,0,0,0,0,0,0,0,0,2,1,0,0,1,0,0,2,0,0,0,1,3,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,2,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,3,0,0,1,1,3,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
1,2,3,1,0,0,2,1,2,1,0,3,0,0,0,0,0,1,1,1,0,1,12,3,5,0,1,3,1,2,5,1,0,1,1,2,0,1,4,0,0,3,3,5,0,0,2,0,0,2,0,0,0,1,4,1,1,0,0,0,1,0,0,0,0,0,2,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,2,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,3,3,4,0,1,1,1,1,3,1,0,0,0,0,0,0,1,0,0,0,2,1,0,0,1,0,0,2,0,0,0,1,3,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,2,0,0,0,0,0,1,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,3,0,0,1,1,3,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
2,6,5,6,0,0,3,2,4,3,5,3,1,0,0,0,0,2,2,2,0,3,16,7,7,3,3,6,4,2,8,4,0,7,2,4,0,2,8,3,2,7,5,7,0,0,3,1,2,2,0,0,0,2,6,4,1,0,1,3,1,0,2,0,0,1,2,3,0,0,1,2,2,1,3,1,1,0,0,0,0,0,0,0,0,0,2,1,0,0,0,1,1,1,3,0,0,2,1,0,0,1,3,1,0,2,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,2,0,0,2,0,0,0,0,0,0,0,0,1,1,0,2,5,4,5,1,2,3,1,1,4,2,0,1,1,1,0,0,1,1,1,2,2,2,0,0,2,0,2,2,0,0,0,2,4,2,0,0,1,2,1,0,0,0,0,2,3,2,0,0,0,0,2,0,2,2,0,0,0,0,0,2,0,0,0,0,5,2,1,1,0,1,1,0,0,2,0,4,0,3,0,0,4,0,1,1,3,3,0,0,0,0,0,0,0,0,0,0,2,0,1,0,0,1,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,4,0,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,1,0,2,0,1,0,0,1,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0
3,6,6,4,1,12,2,0,5,3,9,2,3,1,2,5,0,5,12,8,4,4,12,2,6,0,2,2,0,2,8,8,2,4,0,4,3,10,7,2,5,4,3,3,3,1,0,4,2,2,1,1,0,3,3,0,4,4,10,6,1,0,1,3,2,1,4,0,1,4,1,0,1,0,3,0,1,1,0,0,0,1,4,5,2,0,1,0,1,0,0,0,0,0,2,2,1,0,0,1,2,4,4,0,1,1,0,1,0,0,0,1,2,1,0,0,0,1,1,0,0,0,4,0,1,0,0,0,0,2,1,2,0,8,0,0,4,3,3,1,1,0,1,3,0,1,3,0,0,2,6,1,3,0,0,2,0,1,4,1,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,2,1,0,4,1,3,4,0,0,0,1,1,2,1,1,0,0,0,0,0,0,3,1,0,0,1,2,0,2,5,1,2,2,3,0,2,0,2,0,0,0,1,3,1,2,0,1,1,0,1,1,1,1,1,1,0,0,0,3,0,1,0,0,0,0,0,0,0,2,2,0,0,0,0,2,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,2,0,0,2,1,0,0,0,0,0,1,1,2,0,2,0,2,0,5,2,1,2,2,1,1,2,1,0,0,0,0,0,0,0,0,1,0,0,0,1,2,0,0,1,0,1
4,19,9,1,7,7,8,0,8,4,1,1,3,2,6,5,6,9,8,4,6,4,4,0,8,0,1,2,0,3,8,6,10,3,2,2,1,7,0,3,3,2,2,0,5,4,1,7,3,4,4,2,5,9,5,0,8,2,9,6,7,6,7,6,18,12,4,1,2,1,5,0,3,0,0,0,1,0,4,1,2,2,0,1,4,1,3,0,3,0,0,0,0,0,2,2,4,2,2,0,0,0,0,0,1,0,0,0,1,1,0,2,0,1,1,0,0,2,1,0,3,1,0,2,2,2,1,2,1,4,2,0,3,5,1,0,2,2,1,1,0,0,1,3,2,1,2,1,0,1,0,0,3,0,0,0,0,2,2,1,3,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,1,0,1,1,3,2,0,3,0,3,2,1,2,2,2,4,0,1,0,2,1,0,0,1,1,0,0,1,0,0,0,1,3,3,0,1,0,0,0,0,0,1,2,0,0,1,1,1,1,0,1,0,1,0,1,0,0,0,0,0,1,0,3,1,0,1,0,3,1,0,0,0,0,5,1,1,1,0,0,2,3,1,0,0,0,2,0,2,1,0,0,1,2,1,1,1,3,3,2,1,2,1,0,2,0,0,0,0,1,3,2,2,0,0,1,1,6,0,2,1,2,1,0,3,2,1,2,2,2,2,1,1,3,2,0,2,1,1,1,3,1,4,2,11
'''
    true = pd.read_csv(six.StringIO(s), index_col=0)

    pdt.assert_frame_equal(test, true)


def test_score_kmers(pwm):
    import kvector

    kmers = kvector.make_kmers(3)

    test = pd.Series(kvector.score_kmers(pwm, kmers),
                     index=kmers)
    s = '''AAA,0.4421139794502956
AAC,0.3010679195832866
AAG,0.32337240282664487
AAT,0.42448485149540616
ACA,0.3128897488745921
ACC,0.1718436890075831
ACG,0.19414817225094133
ACT,0.2952606209197027
AGA,0.3128897488745921
AGC,0.1718436890075831
AGG,0.19414817225094133
AGT,0.2952606209197027
ATA,0.5067260947381473
ATC,0.3656800348711383
ATG,0.38798451811449663
ATT,0.48909696678325787
CAA,0.2935687551893772
CAC,0.15252269532236817
CAG,0.17482717856572644
CAT,0.2759396272344877
CCA,0.16434452461367371
CCC,0.023298464746664645
CCG,0.045602947990022916
CCT,0.1467153966587842
CGA,0.16434452461367371
CGC,0.023298464746664645
CGG,0.045602947990022916
CGT,0.1467153966587842
CTA,0.3581808704772289
CTC,0.2171348106102199
CTG,0.23943929385357818
CTT,0.3405517425223395
GAA,0.2935687551893772
GAC,0.15252269532236817
GAG,0.17482717856572644
GAT,0.2759396272344877
GCA,0.16434452461367371
GCC,0.023298464746664645
GCG,0.045602947990022916
GCT,0.1467153966587842
GGA,0.16434452461367371
GGC,0.023298464746664645
GGG,0.045602947990022916
GGT,0.1467153966587842
GTA,0.3581808704772289
GTC,0.2171348106102199
GTG,0.23943929385357818
GTT,0.3405517425223395
TAA,0.4420016217023659
TAC,0.30095556183535693
TAG,0.3232600450787152
TAT,0.4243724937474765
TCA,0.3127773911266625
TCC,0.1717313312596534
TCG,0.1940358145030117
TCT,0.29514826317177295
TGA,0.3127773911266625
TGC,0.1717313312596534
TGG,0.1940358145030117
TGT,0.29514826317177295
TTA,0.5066137369902177
TTC,0.3655676771232087
TTG,0.387872160366567
TTT,0.48898460903532825
'''
    true = pd.read_csv(six.StringIO(s), index_col=0, squeeze=True, header=None)
    true.name = None
    true.index.name = None

    pdt.assert_series_equal(test, true)


def test_score_kmers_kmer_same_length_as_motif(pwm):
    import kvector

    k = 7
    kmers = ['A' * k, 'C' * k, 'G' * k, 'T' * k]

    test = pd.Series(
        kvector.score_kmers(pwm, kmers), index=kmers)

    s = '''AAAAAAA,0.3847637471768403
CCCCCCC,0.04112006099350329
GGGGGGG,0.088915382229271
TTTTTTT,0.48520081057333886
'''
    true = pd.read_csv(six.StringIO(s), index_col=0, squeeze=True, header=None)
    true.index.name = None
    true.name = None

    pdt.assert_series_equal(test, true)


def test_score_kmers_kmer_longer_than_motif(pwm):
    import kvector
    k = 10

    kmers = ['A' * k, 'C' * k, 'G' * k, 'T' * k]

    test = pd.Series(
        kvector.score_kmers(pwm, kmers), index=kmers)

    s = '''AAAAAAAAAA,0.3847637471768403
CCCCCCCCCC,0.04112006099350329
GGGGGGGGGG,0.088915382229271
TTTTTTTTTT,0.48520081057333886
'''
    true = pd.read_csv(six.StringIO(s), index_col=0, squeeze=True, header=None)
    true.index.name = None
    true.name = None

    pdt.assert_series_equal(test, true)


def test_per_interval_kmers(intervals_bed, genome_fasta,
                            interval_kmers_csv):
    import kvector

    test = kvector.per_interval_kmers(intervals_bed, genome_fasta)

    true = pd.read_csv(interval_kmers_csv, index_col=0)

    pdt.assert_frame_equal(test, true)


def test_per_interval_kmers_unnamed_intervals(other_bed, genome_fasta,
                                              other_kmers_csv):
    import kvector

    test = kvector.per_interval_kmers(other_bed, genome_fasta)

    true = pd.read_csv(other_kmers_csv, index_col=0)

    pdt.assert_frame_equal(test, true)


def test_per_interval_kmers_intersect_serial(
        intervals_bed, genome_fasta, other_bed, intervals_intersect_other_kmers_csv):
    import kvector

    test = kvector.per_interval_kmers(intervals_bed, genome_fasta,
                                      intersect=other_bed, threads=0)

    true = pd.read_csv(intervals_intersect_other_kmers_csv, index_col=0)

    pdt.assert_frame_equal(test, true)

def test_per_interval_kmers_intersect_parallelized(
        intervals_bed, genome_fasta, other_bed, intervals_intersect_other_kmers_csv):
    import kvector

    test = kvector.per_interval_kmers(intervals_bed, genome_fasta,
                                      intersect=other_bed, threads=-1)

    true = pd.read_csv(intervals_intersect_other_kmers_csv, index_col=0)

    pdt.assert_frame_equal(test, true)
