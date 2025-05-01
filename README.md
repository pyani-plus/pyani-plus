# `pyANI-plus`

`pyani-plus` is an application and Python module for whole-genome classification of microbes using Average Nucleotide Identity and similar methods. It is a reimplemented version of [`pyani`](https://github.com/widdowquinn/pyani) with support for additional schedulers and methods.

![Linux build](https://github.com/pyani-plus/pyani-plus/actions/workflows/build-linux.yaml/badge.svg)
![macOS build](https://github.com/pyani-plus/pyani-plus/actions/workflows/build-macos.yaml/badge.svg)

[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6b681f069a0443f7b2d7774dbb55de3d)](https://app.codacy.com/gh/pyani-plus/pyani-plus/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/6b681f069a0443f7b2d7774dbb55de3d)](https://app.codacy.com/gh/pyani-plus/pyani-plus/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_coverage)

[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pyani-plus/pyani-plus/main.svg)](https://results.pre-commit.ci/latest/github/pyani-plus/pyani-plus/main)
[![CodeFactor](https://www.codefactor.io/repository/github/pyani-plus/pyani-plus/badge)](https://www.codefactor.io/repository/github/pyani-plus/pyani-plus)
[![codecov](https://codecov.io/gh/pyani-plus/pyani-plus/graph/badge.svg?token=NSSTP6CIW0)](https://codecov.io/gh/pyani-plus/pyani-plus)

## Citing `pyANI-plus`

A complete guide to citing `pyani` is included in the file [`CITATIONS`](CITATIONS). Please cite the following manuscript in your work, if you have found `pyani` useful:

> Pritchard *et al.* (2016) "Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens" *Anal. Methods* **8**, 12-24
DOI: [10.1039/C5AY02550H](https://doi.org/10.1039/C5AY02550H)

## Installation

There are currently no stable releases of `pyani-plus`. If you would like to use the in-progress development version, please follow the usual installation procedure for GitHub repositories, e.g.

1. Clone the repository: `git clone git@github.com:pyani-plus/pyani-plus.git`
2. Change directory to the repository: `cd pyani-plus`
3. Create a new conda environment called `pyani-plus_py312` using the command `make setup_conda_env` (there is a corresponding `remove_conda_env` target)
4. Activate the conda environment with the command `conda activate pyani-plus_py312`
5. Install using one of the following methods:
   1.  `pip`, e.g.: `pip install -U -e .`
   2.  `Make`, e.g.: `make install_macos` or `make install_linux`

## Walkthrough: A First Analysis

## Contributing

Please see the [`CONTRIBUTING.md`](CONTRIBUTING.md) file for more information

## Method and Output Description

## Licensing

Unless otherwise indicated, the material in this project is made available under the MIT License.

```text
    (c) The University of Strathclyde 2024-present
    Contact: leighton.pritchard@strath.ac.uk

    Address:
    Dr Leighton Pritchard,
    Strathclyde Institute of Pharmacy and Biomedical Sciences
    161 Cathedral Street
    Glasgow
    G4 0RE,
    Scotland,
    UK

The MIT License

Copyright (c) 2024-present The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```
