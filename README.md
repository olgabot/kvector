# kvector

[![](https://img.shields.io/travis/olgabot/kvector.svg)](https://travis-ci.org/olgabot/kvector)[![](https://img.shields.io/pypi/v/kvector.svg)](https://pypi.python.org/pypi/kvector)

## What is `kvector`?

kvector is a small utility for converting motifs to kmer vectors to compare motifs of different lengths

* Free software: BSD license
* Documentation: https://olgabot.github.io/kvector

## Installation

To install this code, clone this github repository and use `pip` to install

    git clone git@github.com:olgabot/kvector
    cd kvector
    pip install .  # The "." means "install *this*, the folder where I am now"


## Features

### Read HOMER motif file

```python
motifs = kvector.read_motifs("kvector/tests/example.motifs", residues='ACGT')
```

The output is a pandas Series of the motif ids from the file, mapped to a dataframe of the position-weight matrix of the motif.

### Create metadata matrix from the ID lines of the motifs

```python
metadata = kvector.create_metadata(motifs)
```

### Transform the motif PWM to a kmer vector

```python
motif_kmer_vectors = kvector.motifs_to_kmer_vectors(motifs, residues='ACGT', 
    kmer_lengths=kmer_lengths)
```