import os
import pytest


@pytest.fixture(params=((3, 4), 3))
def kmer_lengths(request):
    return request.param


@pytest.fixture
def data_folder():
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


@pytest.fixture
def intervals_bed(data_folder):
    return os.path.join(data_folder, 'intervals.bed')


@pytest.fixture
def other_bed(data_folder):
    return os.path.join(data_folder, 'other.bed')


@pytest.fixture
def genome_fasta(data_folder):
    return os.path.join(data_folder, 'chromosome.fasta')
