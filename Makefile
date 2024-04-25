# Makefile
#
# This file is part of the pyani-plus package distribution
# https://github.com/pyani-plus/pyani-plus

# Install OS-independent conda dependencies
setup_conda:
	@conda install --file requirements-dev.txt --yes

# Set up development environment (OS-dependent)
setup_dev_linux: setup_conda
	@pre-commit install
	@pip install -U -e .

setup_dev_macos: setup_conda
	@pre-commit install
	@pip install -U -e .

# Install pyani-plus (OS-dependent)
install_linux: setup_conda
	@pip install -U -e .

install_macos: setup_conda
	@pip install -U -e .

# Run tests
# When the tests complete, the coverage output will be opened in a browser
test:
	@python -m pytest --cov-report=html --cov=pyani-plus -v tests/ && open htmlcov/index.html
