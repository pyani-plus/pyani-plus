# Makefile
#
# This file is part of the pyani-plus package distribution
# https://github.com/pyani-plus/pyani-plus

# Install OS-independent conda dependencies
setup_conda:
	@conda install --file requirements.txt --yes

# Install OS-independent conda dependencies for development
setup_conda-dev:
	@conda install --file requirements-dev.txt --yes

# Install pyani-plus (OS-dependent, not developer version)
install_linux: setup_conda
	@conda install --file requirements-thirdparty-linux.txt --yes
	@pip install -U -e .

install_macos: setup_conda
	@conda install --file requirements-thirdparty-macos.txt --yes
	@pip install -U -e .

# Set up development environment (OS-dependent)
setup_dev_linux: setup_conda setup_conda-dev
	@conda install --file requirements-thirdparty-linux.txt --yes
	@pre-commit install
	@pip install -U -e .

setup_dev_macos: setup_conda setup_conda-dev
	@conda install --file requirements-thirdparty-macos.txt --yes
	@pre-commit install
	@pip install -U -e .

# Run tests
# When the tests complete, the coverage output will be opened in a browser
test:
	@python -m pytest -n auto --cov-report=html --cov=pyani_plus -v tests/ && open htmlcov/index.html

# Clean up test output
clean_test:
	@rm -rf htmlcov
	@rm -rf tests/nucmer_filter_output
	@rm -rf tests/nucmer_delta_output
