import six
import pandas as pd

from .kmer import score_kmers, make_kmers, RNA, DNA


def homer_motif_reader(handle, residues=DNA):
    """Read homer motifs output and return tuples of motif_id, motif_pwm

    Parameters
    ----------
    handle : file
        A Python opened file object
    residues : str
        Name of the residues
    """
    names = list(residues)
    record_id, record = None, ''
    for line in handle:
        if line.startswith('>'):
            new_record_id = line.lstrip('>').strip()
            if record_id is None:
                record_id = new_record_id
            if len(record) > 0:
                pwm = pd.read_table(six.StringIO(record), header=None, names=names)
                yield record_id, pwm
                record = ''
                record_id = new_record_id
        else:
            record += line



def read_motifs(filename, residues):
    """Wrapper to read a homer motif file

    Parameters
    ----------
    filename : str
        Name of the motif file to open
    residues : str

    Returns
    -------
    motifs : pandas.Series
        A series of dicts holding the name of the motif

    """
    with open(filename) as f:
        motifs = pd.Series(dict(homer_motif_reader(f, residues=residues)))
    return motifs

def create_metadata(motifs):
    motif_metadata = pd.DataFrame.from_records(
        motifs.index.map(lambda x: pd.Series(x.split())))

    motif_metadata.columns = ['Consensus Sequence', 'Motif ID',
                              'Log Odds Threshold', '$\log(p)$',
                              'empty_placeholder', 'Occurence Information',
                              'Motif Statistics']
    return motif_metadata

def motifs_to_kmer_vectors(motifs, residues):
    kmer_lengths = motifs['Consensus Sequence'].map(len).unique()
    kmers = make_kmers(kmer_lengths, residues)
    kmers_list = map(list, kmers)

    motif_scores = motifs.map(
        lambda x: pd.Series(score_kmers(x, kmers_list), index=kmers))
    motif_scores = pd.DataFrame.from_records(motif_scores).T
    motif_scores.columns = motifs.index
    return motif_scores
