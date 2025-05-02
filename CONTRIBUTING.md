# Contributing to `pyani-plus`

Contributions, including bugfixes and addition of new features, are welcomed!

This document provides a brief introduction to contributing to `pyani-plus`. More complete documentation will be made available elsewhere.

## Licensing

`pyani-plus` is licensed under a permissive MIT Licence (see `LICENSE` file for details). Any contributions you make will be licensed under this agreement.

## Repository branches

The `pyani-plus` package is maintained on [GitHub](https://github.com/pyani-plus/pyani-plus). The current development version is always maintained under the `main` branch. In addition, we keep the following "live" branches for convenience:

- `github_administration`: intended for development/modification of GitHub-related files (e.g. GitHub Actions, templates, etc.)
- `dev_administration`: intended for development/modification of developer-related activities (e.g. configuration of linters, testing, etc.)

We also develop in dynamic branches, which may be created/destroyed over time. For instance, if working on issue #49, we would manage this on a branch named `issue_49` and supporting scratch work should take place in the subdirectory `issue_49`. Note that all subdirectories beginning `issue_*` are ignored by a rule in `.gitignore`.

Other branches may be created as and when necessary. Please name branches clearly and unambiguously.

## How to Contribute

### Clone the repository

If you are not in the `pyani-plus` development team, please fork this repository to your own GitHub account (you will need to create a GitHub account if you do not already have one), and clone the repository from your own fork:

```bash
git clone git@github.com:<YOUR USERNAME>/pyani-plus.git
cd pyani-plus
```

### Set up the development environment

We recommend creating a `conda` environment specifically for development of `pyani-plus`, and use just the bioconda and conda-forge channels:

```bash
conda create --name pyani-plus_py312 python=3.12 -y
conda activate pyani-plus_py312
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
conda config --remove channels defaults
```

With the environment activated, use the `Makefile` to set up the recommended development environment. Use the `Make` rule corresponding to your operating system

```bash
make setup_dev_linux  # Linux OR...
make setup_dev_macos  # macOS
```

Those conda and make commands are a one-off setup, after which all you need in a fresh terminal session is:

```bash
conda activate pyani-plus_py313
```

This will set up the tool pre-commit as a git pre-commit hook, which will run
assorted checks including ruff for formatting and style checking.

This will also install development dependencies like ``pytest`` which we use to run
our test suite. Run the test with ``make test``, or by directly calling pytest if
you want to modify the arguments to pytest, e.g. ``pytest -v -n auto`` will run with
multiple worker threads which can be significantly faster on a many-core machine.

### Commit message conventions

`git` commit messages are an important way to manage a readable revision history. We use the following conventions:

- If a commit fixes an issue, state this in the commit message
  - GitHub Projects picks up on terms like `fixes #123` and `closes #987`. Using these phrases makes project management much easier.

**Every commit gets a short description**

- The short description should be in imperative form and around 50 characters or less
- The short description can, but need not, include the name of the file that is modified
- There should be a short account of what the change does
- The short description should not only contain the name of the file that is modified
- The short description should continue the sentence "This commit will..."

For example, if the commit updates some documentation, the following are good short descriptions:

- `update citations.rst to add new references`
- `update docs to add new references; fixes #567`
- `add new references to citations.rst`

The following are not good short descriptions

- `update citations.rst` (does not say what was done)
- `there were some new references so I added them` (not in imperative form)
- `citations.rst` (does not say what was done)
- `part of some doc updates` (does not say what was done)

**Some commits get long/extended descriptions**

- Long descriptions should be included where there is more information to share than can be fit in the short description
- Long descriptions are free-form
- Long descriptions should explain the why of the change. They may also explain the what at a high level, but should not be excessively detailed.
- They are "long descriptions", but they should be concise, precise and clear
- Paragraphs are fine
- Bullet points are fine

### Contributing your changes

Please use a short, descriptive branch name to make your changes in. If you are addressing an issue from the [Issues](https://github.com/pyani-plus/pyani-plus/issues) page, please use the branch name `issue_N` where `N` is the number of the issue. Once you have finished making your changes, please push them to your fork and submit a pull request against the `pyani-plus` repository.

A typical workflow might look something like this:

1. Identify something you think you can change or add
2. Fork this repository into your own account
3. Obtain source code and set up a development environment as outlined above
4. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
5. When you're finished, test locally with `pytest -v` or `make test`.
6. Push the branch to your fork, and submit a pull request (please tick "Allow edits by maintainers")
7. Continue the discussion in the [Pull Requests](https://github.com/pyani-plus/pyani-plus/pulls) section of this repository on GitHub.

Release process
---------------

For a release, start from a clean git checkout. You will need some python tools not covered by the above development setup:

```bash
pip install -U pip twine build
```

After checking the version defined in `pyani_plus/__init__.py` and `README.md` are up to date (and checked in):

```bash
rm -rf build/
python -m build
git tag vX.Y.Z
git push origin main --tags
twine upload dist/pyani_plus-X.Y.Z*
```

Then you must also turn the git tag into a "release" on GitHub: https://github.com/pyani-plus/pyani-plus/releases
This should automatically generate a version specific DOI on Zenodo: https://doi.org/10.5281/zenodo.15005805
