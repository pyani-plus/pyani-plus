# `pyANI-plus`

`pyani-plus` is an application and Python module for whole-genome classification of microbes using Average Nucleotide Identity and similar methods. It is a reimplemented version of [`pyani`](https://github.com/widdowquinn/pyani) with support for additional schedulers and methods.

![Linux build](https://github.com/pyani-plus/pyani-plus/actions/workflows/build-linux.yaml/badge.svg)
![macOS build](https://github.com/pyani-plus/pyani-plus/actions/workflows/build-macos.yaml/badge.svg)

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6b681f069a0443f7b2d7774dbb55de3d)](https://app.codacy.com/gh/pyani-plus/pyani-plus/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/6b681f069a0443f7b2d7774dbb55de3d)](https://app.codacy.com/gh/pyani-plus/pyani-plus/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pyani-plus/pyani-plus/main.svg)](https://results.pre-commit.ci/latest/github/pyani-plus/pyani-plus/main)
[![CodeFactor](https://www.codefactor.io/repository/github/pyani-plus/pyani-plus/badge)](https://www.codefactor.io/repository/github/pyani-plus/pyani-plus)
[![codecov](https://codecov.io/gh/pyani-plus/pyani-plus/graph/badge.svg?token=NSSTP6CIW0)](https://codecov.io/gh/pyani-plus/pyani-plus)

[![PyPI](https://img.shields.io/pypi/v/pyani_plus.svg?label=PyPI)](https://pypi.org/project/pyani-plus/)
[![BioConda](https://img.shields.io/conda/vn/bioconda/pyani-plus.svg?label=Bioconda)](https://anaconda.org/bioconda/pyani-plus)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15005805.svg)](https://doi.org/10.5281/zenodo.15005805)
[![Code Style](https://img.shields.io/badge/Code%20style-black-000000.svg)](https://github.com/python/black)

## Citing `pyANI-plus`

A complete guide to citing `pyani` is included in the file [`CITATIONS`](CITATIONS). Please cite the following manuscript in your work, if you have found `pyani` useful:

> Pritchard *et al.* (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" *Anal. Methods* **8**, 12-24
DOI: [10.1039/C5AY02550H](https://doi.org/10.1039/C5AY02550H)

## Documentation

See [pyANI-plus documentation](https://pyani-plus.github.io/pyani-plus-docs/) including [Walkthrough: A First Analysis](https://pyani-plus.github.io/pyani-plus-docs/walkthrough.html).

## Installation

See our [full installation instructions](https://pyani-plus.github.io/pyani-plus-docs/installation.html) for more details. In brief, if you are using `conda` we recommend
[`pyani-plus` from the BioConda channel](https://anaconda.org/bioconda/pyani-plus),
which will include all the dependencies. Otherwise, we recommend installing [`pyani-plus`
from the Python Package Index (PyPI)](https://pypi.org/project/pyani-plus/),
and install the non-python tools as needed.

If you would like to use the in-progress development version, please follow the usual installation procedure for GitHub repositories, e.g.

1. Clone the repository: `git clone git@github.com:pyani-plus/pyani-plus.git`
2. Change directory to the repository: `cd pyani-plus`
3. Create a new conda environment called `pyani-plus_py312` using the command `make setup_conda_env` (there is a corresponding `remove_conda_env` target)
4. Activate the conda environment with the command `conda activate pyani-plus_py312`
5. Install using one of the following methods:
   1.  `pip`, e.g.: `pip install -U -e .`
   2.  `Make`, e.g.: `make install_macos` or `make install_linux`

## Contributing

Please see the [`CONTRIBUTING.md`](CONTRIBUTING.md) file for more information,
including installing additional development-only dependencies like `pytest` and `ruff`.

## Licensing

Unless otherwise indicated, the material in this project is made available under the MIT License,
and copyright The University of Strathclyde 2024 to present. See the LICENSE file.

> Contact: leighton.pritchard@strath.ac.uk
>
> Address:
>    Dr Leighton Pritchard,
>    Strathclyde Institute of Pharmacy and Biomedical Sciences
>    161 Cathedral Street
>    Glasgow
>    G4 0RE,
>    Scotland,
>    UK
