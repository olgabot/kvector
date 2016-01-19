#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_kmer
----------------------------------

Tests for `kvector.kmer` module.
"""

import pandas as pd
import pandas.util.testing as pdt
import pytest
import six


@pytest.fixture
def fasta(data_folder):
    return '{}/example.fasta'.format(data_folder)

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
