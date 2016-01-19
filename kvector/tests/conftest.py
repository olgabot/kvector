import os
import pytest


@pytest.fixture(params=((3, 4), 3))
def kmer_lengths(request):
    return request.param


@pytest.fixture
def data_folder():
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
