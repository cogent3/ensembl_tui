[![CI](https://github.com/cogent3/ensembl_tui/actions/workflows/testing_develop.yml/badge.svg)](https://github.com/cogent3/ensembl_tui/actions/workflows/testing_develop.yml)
[![CodeQL](https://github.com/cogent3/ensembl_tui/actions/workflows/codeql.yml/badge.svg)](https://github.com/cogent3/ensembl_tui/actions/workflows/codeql.yml)
[![Coverage Status](https://coveralls.io/repos/github/cogent3/ensembl_tui/badge.svg?branch=develop)](https://coveralls.io/github/cogent3/ensembl_tui?branch=develop)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![DOI](https://zenodo.org/badge/664119735.svg)](https://doi.org/10.5281/zenodo.15098645)
[![DOI](https://app.readthedocs.org/projects/ensembl-tui/badge/?version=latest)](https://ensembl-tui.readthedocs.io/en/latest/)

# ensembl-tui

ensembl-tui provides the `eti` terminal application for obtaining a subset of the data provided by Ensembl which can then be queried locally. You can have multiple such subsets on your machine, each corresponding to a different selection of species and data types.

> **Warning**
> We currently **only support accessing data from the main ensembl.org** site. If you discover errors, please post a [bug report](https://github.com/cogent3/ensembl_tui/issues).

## Installing the software

<details>
  <summary>General user installation instructions</summary>

  ```
  $ pip install ensembl-tui
  ```

</details>

<details>
  <summary>Developer installation instructions</summary>
  Fork the repo and clone your fork to your local machine. In the terminal, create either a python virtual environment or a new conda environment and activate it. In that virtual environment

  ```
  $ pip install flit
  ```

  Then do the flit version of a "developer install". (It is basically creating a symlink to the repos source directory.)

  ```
  $ flit install -s --python `which python`
  ```
</details>

## Resources required to subset Ensembl data

Ensembl hosts some very large data sets. You need to have a machine with sufficient disk space to store the data you want to download. At present we do not have support for predicting how much storage would be required for a given selection of species and data types. You will need to experiment.

Some commands can be run in parallel but have moderate memory requirements. If you have a machine with limited RAM, you may need to reduce the number of parallel processes. Again, run some experiments.

## Getting setup

<details>
  <summary>Specifying what data you want to download and where to put it</summary>

  We use a plain text file to indicate the Ensembl domain, release and types of genomic data to download. Start by using the `demo-config` subcommand.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["demo-config", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti demo-config [OPTIONS]

    exports sample config and species table to the nominated path

  Options:
    -o, --outpath PATH     Path to directory to export all rc contents.
    -f, --force_overwrite  Overwrite existing data.
    --help                 Show this message and exit.

  ```
  <!-- [[[end]]] -->

  ```shell
  $ eti demo-config -o ensembl_download
  ```
  This command creates a `ensembl_download` download directory and writes two plain text files into it:

  1. `species.tsv`: contains the Latin names, common names etc... of the species accessible at ensembl.org website.
  2. `sample.cfg`: a sample configuration file that you can edit to specify the data you want to download.

  The latter file includes comments on how to edit it in order to specify the genomic resources that you want.

</details>

<details>
  <summary>Downloading the data</summary>

  Downloads the data indicated in the config file to a local directory.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["download", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti download [OPTIONS]

    download data from Ensembl's ftp site

  Options:
    -c, --configpath PATH  Path to config file specifying databases, (only species
                           or compara at present).
    -d, --debug            Maximum verbosity, and reduces number of downloads,
                           etc...
    -v, --verbose
    --help                 Show this message and exit.

  ```
  <!-- [[[end]]] -->

  For a config file named `config.cfg`, the download command would be:

  ```shell
  $ cd to/directory/with/config.cfg
  $ eti download -c config.cfg
  ```

  > **Note**
  > This is the only step for which the internet is required. Downloads can be interrupted and resumed. The software will delete partially downloaded files.

The download creates a new `.cfg` file inside the download directory. This file is used by the `install` command.

</details>

<details>
  <summary>Installing the data</summary>

Converts the downloaded data into data formats designed to enhance querying performance.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["install", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti install [OPTIONS]

    create the local representations of the data

  Options:
    -d, --download PATH       Path to local download directory containing a cfg
                              file.
    -np, --num_procs INTEGER  Number of procs to use.  [default: 1]
    -f, --force_overwrite     Overwrite existing data.
    -v, --verbose
    --help                    Show this message and exit.

  ```
  <!-- [[[end]]] -->

This step can be run in parallel, but the memory requirements will scale with the number of genomes. So we suggest monitoring performance on your system by trying it out on a small number of CPUs to start with. The following command uses 2 CPUs and has been safe on systems with only 16GB of RAM for 10 primate genomes, including homology data and whole genome alignments.

```shell
$ cd to/directory/with/downloaded_data
$ eti install -d downloaded_data -np 2
```

</details>

<details>
  <summary>Checking what has been installed</summary>
This will give a summary of what data has been installed at a provided path.


  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main, ["installed", "--help"])
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti installed [OPTIONS]

    show what is installed

  Options:
    -i, --installed TEXT  Path to root directory of an installation.  [required]
    --help                Show this message and exit.

  ```
  <!-- [[[end]]] -->

</details>

## Interrogating the data

We provide a conventional command line interface for querying the data with subcommands.

<details>
  <summary>The full list of subcommands</summary>

  You can get help on individual subcommands by running `eti <subcommand>` in the terminal.

  <!-- [[[cog
  import cog
  from ensembl_tui import cli
  from click.testing import CliRunner
  runner = CliRunner()
  result = runner.invoke(cli.main)
  help = result.output.replace("Usage: main", "Usage: eti")
  cog.out(
      "```\n{}\n```".format(help)
  )
  ]]] -->
  ```
  Usage: eti [OPTIONS] COMMAND [ARGS]...

    Tools for obtaining and interrogating subsets of https://ensembl.org genomic
    data.

  Options:
    --version  Show the version and exit.
    --help     Show this message and exit.

  Commands:
    tui              Open Textual TUI.
    demo-config      exports sample config and species table to the nominated...
    download         download data from Ensembl's ftp site
    install          create the local representations of the data
    installed        show what is installed
    species-summary  genome summary data for a species
    dump-genes       export meta-data table for genes from one species to...
    compara-summary  summary data for compara
    homologs         exports CDS sequence data in fasta format for homology...
    alignments       export multiple alignments in fasta format for named genes

  ```
  <!-- [[[end]]] -->

</details>

We also provide an experiment terminal user interface (TUI) that allows you to explore the data in a more interactive way. This is invoked with the `tui` subcommand.

### Getting a summary of a genome

A command like the following
```
eti species-summary -i primates10_113/install --species human
```
displays two tables for the indicated genome. The first is the biotypes and their counts, the second the repeat classes / types and their counts.

### Getting a summary of a homology data

A command like the following
```
eti compara-summary -i primates10_113/install
```
displays the homology types and counts. The values under `homology_type` can be used as input arguments to the `homologs` command `--homology_type` argument.

### Exporting related sequences

A command like the following
```
eti homologs -i primates10_113/install/ --outdir sampled_100 --ref human --coord_names 1 --limit 100
```
will sample 100 one-to-one orthologs (the default homology type) to human chromosome 1 linked protein coding genes (the only biotype supported at present). The canonical CDS sequences will be written in fasta format to the directory `sampled_100`.

### Exporting whole genome alignments

A command like the following
```
eti alignments -i primates10_113/install --outdir sampled_aligns_100 --align_name '*primate*' --coord_names 1 --ref human --limit 10
```
samples 10 alignments that include human chromosome 1 protein coding genes. These are from the Ensembl whole genome alignment whose name matches the glob pattern `*primate*`.

> **Warning**
>
> If this pattern matches more than one installed Ensembl alignment, the program will exit.