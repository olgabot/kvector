kvector
=======

|image0|\ |image1|

What is ``kvector``?
--------------------

kvector is a small utility for converting motifs to kmer vectors to
compare motifs of different lengths

-  Free software: BSD license
-  Documentation: https://olgabot.github.io/kvector

Installation
------------

To install this code, clone this github repository and use ``pip`` to
install

::

    git clone https://github.com/olgabot/kvector.git
    cd kvector
    pip install .  # The "." means "install *this*, the folder where I am now"

Features
--------

| Check out `this
  notebook <https://github.com/olgabot/kvector/blob/master/overview.ipynb>`__
| for an overview of features with both inputs and outputs (below shows
  only inputs)

Count k-mers for each line in a ``bed`` file (multithreaded)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For each interval in a bed file, count the kmers and return a
| (n\_intervals, n\_kmers) matrix of the k-mer counts of each region.

.. code:: python

    kmers = kvector.per_interval_kmers(bedfile, genome_fasta, threads=threads,
        kmer_lengths=(4, 5, 6), residues='ACGT')
    csv = bedfile.replace('.bed', '_kmers.csv')
    kmers.to_csv(csv)

Count k-mers for each line in a ``bed`` file, intersected (multithreaded)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For each interval in a bed file, intersect it with another (``other``)
  bed file (e.g. only
| conserved regions of introns) and count k-mers for the intersected
  region. Returns
| a (n\_intervals, n\_kmers) matrix of the k-mer counts of each line in
  the bed file,
| intersected with the ``other`` bed.

.. code:: python

    kmers = kvector.per_interval_kmers(bedfile, genome_fasta, other, threads=threads,
        kmer_lengths=(4, 5, 6))
    csv = bedfile.replace('.bed', '_kmers.csv')
    kmers.to_csv(csv)

Count all *k*-mers in a fasta file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    kmer_vector = kvector.count_kmers('kvector/tests/data/example.fasta', kmer_lengths=(3, 4))
    kmer_vector.head()

Read HOMER motif file
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    motifs = kvector.read_motifs("kvector/tests/example.motifs", residues='ACGT')

| The output is a pandas Series of the motif ids from the file, mapped
  to a
| dataframe of the position-weight matrix of the motif.

Create metadata matrix from the ID lines of the motifs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    metadata = kvector.create_metadata(motifs)

Transform the motif PWM to a kmer vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Keep in mind that on most computers, only kmers up to about 8 (4^8 =
  65,536)
| can be stored in memory. You may want to do this on a supercomputer
  and not
| just your laptop.

.. code:: python

    motif_kmer_vectors = kvector.motifs_to_kmer_vectors(motifs, residues='ACGT',
        kmer_lengths=(4, 5, 6))

.. |image0| image:: https://img.shields.io/travis/olgabot/kvector.svg
   :target: https://travis-ci.org/olgabot/kvector
.. |image1| image:: https://img.shields.io/pypi/v/kvector.svg
   :target: https://pypi.python.org/pypi/kvector
