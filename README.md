# On the spectrum of mesons in quenched Sp(2N) gauge theories---analysis workflow release

This repository contains the analysis code used to prepare the plots and tables
included in [On the spectrum of mesons in quenched Sp(2N) gauge theories][quenchedmesons-paper].

## Requirements

The analysis is written in a mix of WolframScript and Python. The various
components are joined together using Make. The code has been tested with

* Mathematica 14.0 with MaTeX installed
* [Snakemake 8.2.1][snakemake], installed via Mamba
* Python 3.12
  (with requirements as documented in `environment.yml`,
  and installed using Snakemake)

## Setup

### Installation

* Download the repository
* Ensure that the `wolframscript` is available in your `PATH`.
* Ensure that MaTeX is installed.

      wolframscript -file extra_tools/install_matex.wls

### Data

* Download and extract the data, [available from Zenodo][data] [TODO].
  Specifically,
  * Download the file `raw_data.zip`
    and unzip its contents into the root directory
    to create the `raw_data` subdirectory.
  * Download the file `metadata.zip`
    and unzip its contents into the `metadata` subdirectory.

### Running the analysis

* Verify how many Mathematica licenses you have access to by running

      wolframscript -file extra_tools/check_licenses.wls

  Close all Mathematica instances,
  or subtract one from this number
  to get the number of WolframScript instances you can start.

* With the software and data downloaded,
  you should be able to reproduce the full analysis
  by typing

      snakemake --use-conda --cores [cores] --resource mathematica_licenses=[licenses] --retries 2

  where:

  * `[cores]` is the number of CPU cores you want to allow Snakemake to use, and
  * `[licenses]` is the number of available Mathematica licenses.

  The `--retries` option is needed
  because sometimes Mathematica will fail for no specific reason.

* On eight Apple M1 CPU cores,
  the full analysis takes around an hour and a half.


## Output

The following classes of output are created:

* The `plots` directory contains all plots included in the paper.
* The `tables` directory contains
  LaTeX source for all numerical tables included in the paper.
* The file `definitions.tex` defines macros for
  all numbers quoted in the paper
  that are outputs of the analysis.
* The `csvs` directory contains CSV files,
  including all data presented in plots or tables,
  for ease of those who wish to build on our results in future work.
  Pre-generated versions of these
  are available from the data release. <!-- TODO ADD LINK -->
* The `processed_data` directory contains
  intermediary files not included in the publication.

All directories are created automatically.


## Reproducibility and reusability

Some portions of this analysis make use of
Mathematica functions that are not bitwise reproducible between CPU architectures.
As such,
not all results obtained will be bitwise identical
when run on different computers.
Such differences are well outside the precision we consider in this work,
and do not affect ourconclusions.

This analysis was originally performed manually using Mathematica,
and was automated after the fact
with the specific dataset linked above in mind.
The code included in this repository is therefore
not directly suited for application to other datasets
without some work to either adopt a very similar schema,
or to make the workflow more generalisable.

[quenchedmesons-paper]: https://arxiv.org/abs/2312.08465
[snakemake]: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
