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

If you are not in the `pyani-plus` development team, please fork this repository to your own GitHub account (you will need to create a GitHub account if you do not alread have one), and clone the repository from your own fork:

```bash
git clone git@github.com:<YOUR USERNAME>/pyani-plus.git
cd pyani-plus
```

### Set up the development environment

We recommend creating a `conda` environment specifically for development of `pyani-plus`.

```bash
conda create --name pyani-plus_py312 python=3.12 -y   
conda activate pyani-plus_py312
```

With the environment activated, use the `Makefile` to set up the recommended development environment. Use the `Make` rule corresponding to your operating system

```bash
make setup_dev_linux  # Linux OR...
make setup_dev_macos  # macOS
```

### Commit message conventions

`git` commit messages are an important way to manage a readable revision history. We use the following conventions:

- If a commit fixes an issue, state this in the commit message
  - GitHub Projects picks up on terms like `fixes #123` and `closes #987`. Using these phrases makes project management much easier.

**Every commit gets a short description**

- The short description should be in imperative form and around 50 characters or less
- The short description can, but need not, include the name of the file that is modified
- There should be a short account of what the change does
- The short description should not only contain the name of the file that is modified
- The short description shoudl continue the sentence "This commit will..."

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
2. Fork this repository into your own account (**please add write access for `pyani-plus` maintainers**)
3. Obtain source code and set up a development environment as outlined above
4. Create a new branch with a short, descriptive name (for the thing you're fixing/adding), and work on this branch locally
5. When you're finished, test locally with `pytest -v` or `make test`.
6. Push the branch to your fork, and submit a pull request
7. Continue the discussion in the [Pull Requests](https://github.com/pyani-plus/pyani-plus/pulls) section of this repository on GitHub.

