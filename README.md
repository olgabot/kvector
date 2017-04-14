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

Check out [this notebook](https://github.com/olgabot/kvector/blob/master/overview.ipynb)
for an overview of features with both inputs and outputs (below shows only inputs)

### Count k-mers for each line in a `bed` file (multithreaded)

For each interval in a bed file, count the kmers and return a
(n_intervals, n_kmers) matrix of the k-mer counts of each region.

```python
kmers = kvector.per_interval_kmers(bedfile, genome_fasta, threads=threads,
    kmer_lengths=(4, 5, 6), residues='ACGT')

# This is a pandas DataFrame so you can look at the top ("head") of it
kmers.head()

# Save the
csv = bedfile.replace('.bed', '_kmers.csv')
kmers.to_csv(csv)
```

### Count k-mers for each line in a `bed` file, intersected (multithreaded)

For each interval in a bed file, intersect it with another (`other`) bed file (e.g. only
conserved regions of introns) and count k-mers for the intersected region. Returns
a (n_intervals, n_kmers) matrix of the k-mer counts of each line in the bed file,
intersected with the `other` bed.

```python
kmers = kvector.per_interval_kmers(bedfile, genome_fasta, other, threads=threads,
    kmer_lengths=(4, 5, 6))
csv = bedfile.replace('.bed', '_kmers.csv')
kmers.to_csv(csv)
```

### Count all *k*-mers in a fasta file

```python
kmer_vector = kvector.count_kmers('kvector/tests/data/example.fasta', kmer_lengths=(3, 4))
kmer_vector.head()
```

### Read HOMER motif file

```python
motifs = kvector.read_motifs("kvector/tests/example.motifs", residues='ACGT')
```

The output is a pandas Series of the motif ids from the file, mapped to a
dataframe of the position-weight matrix of the motif.

### Create metadata matrix from the ID lines of the motifs

```python
metadata = kvector.create_metadata(motifs)
```

### Transform the motif PWM to a kmer vector

Keep in mind that on most computers, only kmers up to about 8 (4^8 = 65,536)
can be stored in memory. You may want to do this on a supercomputer and not
just your laptop.

```python
motif_kmer_vectors = kvector.motifs_to_kmer_vectors(motifs, residues='ACGT',
    kmer_lengths=(4, 5, 6))
```

## Running from the command line

You can run `kvector` from the command line. For example, from this git
repository you could run:

```
kvector kvector/tests/data/intervals.bed kvector/tests/data/chromosome.fasta
```

Here's the full help:

```
$ kvector --help
Usage: kvector [OPTIONS] BED FASTA

  Counts k-mers in the bed intervals and writes a csv to stdout

  Parameters
  ----------
  bed : str
      Location of a bed file of the genomic intervals whose kmers you want
      to count
  fasta : str
      Path to the genome fasta file containing all chromosomes. Must be
      indexed (usually has a `.fai` file in the same directory, created
      using `faidx`).

Options:
  --intersect TEXT     Bed file of other regions, e.g. conserved elements,
                       that you want to intersect when searching for k-mers.
  --kmer-lengths TEXT  How long of DNA words to search for (aka values of k).
                       Default is "4,5,6".
  --residues TEXT      Which letters to search for in the fasta file. Default
                       is 'ACGT'.
  --threads INTEGER    Number of threads/processors to use for paralell
                       processing of a multithreaded job. Default is -1, which
                       uses the maximum number of threads available, via the
                       "joblib" module.
  --version            Show the version and exit.
  --help               Show this message and exit.
```

## How to run the tests

### Vanilla test running

To run the tests and nothing else, type:

```
py.test
```

### Tests with coverage

To calculate how many lines of code are covered by tests, use the `Makefile`
command:

```
make coverage
```
