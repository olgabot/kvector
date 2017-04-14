import subprocess

from click.testing import CliRunner


def test_cli_help():
    from kvector.commandline import cli

    runner = CliRunner()
    test = runner.invoke(cli, ['--help'])

    true_output = """Usage: cli [OPTIONS] BED FASTA

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
  --intersect TEXT     Bed file of other regions, e.g. conserved elements, that
                       you want to intersect when searching for k-mers.
  --kmer-lengths TEXT  How long of DNA words to search for (aka values of k).
                       Default is "4,5,6".
  --residues TEXT      Which letters to search for in the fasta file. Default is
                       'ACGT'.
  --threads INTEGER    Number of threads/processors to use for paralell
                       processing of a multithreaded job. Default is -1, which
                       uses the maximum number of threads available, via the
                       "joblib" module.
  --version            Show the version and exit.
  --help               Show this message and exit.
"""

    assert test.output == true_output


def test_kvector_installed(capsys):
    test_output = subprocess.call(['kvector'])
    assert False


def test_cli_with_arguments(intervals_bed, genome_fasta):
    from kvector.commandline import cli

    runner = CliRunner()
    test = runner.invoke(cli, [intervals_bed, genome_fasta])

    true_output = """Usage: cli [OPTIONS] BED FASTA

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
  --intersect TEXT     Bed file of other regions, e.g. conserved elements, that
                       you want to intersect when searching for k-mers.
  --kmer-lengths TEXT  How long of DNA words to search for (aka values of k).
                       Default is "4,5,6".
  --residues TEXT      Which letters to search for in the fasta file. Default is
                       'ACGT'.
  --threads INTEGER    Number of threads/processors to use for paralell
                       processing of a multithreaded job. Default is -1, which
                       uses the maximum number of threads available, via the
                       "joblib" module.
  --version            Show the version and exit.
  --help               Show this message and exit.
"""

    assert test.output == true_output
