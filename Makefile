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

fixtures:
	@echo "Running scripts to recreate the input files to the test suite."
	@echo "This will take several minutes..."
	@echo ""
	cd tests/generate_fixtures/; ./generate_anib_fragment_files.py
	cd tests/generate_fixtures/; ./generate_anib_blast_files.py
	cd tests/generate_fixtures/; ./generate_target_anim_files.py ../fixtures/viral_example ../fixtures/viral_example/intermediates/ANIm
	cd tests/generate_fixtures/; ./generate_target_anim_files.py ../fixtures/bad_alignments ../fixtures/bad_alignments/intermediates/ANIm
	cd tests/generate_fixtures/; ./generate_target_dnadiff_files.py ../fixtures/viral_example ../fixtures/viral_example/intermediates/dnadiff
	cd tests/generate_fixtures/; ./generate_target_dnadiff_files.py ../fixtures/bad_alignments ../fixtures/bad_alignments/intermediates/dnadiff
	cd tests/generate_fixtures/; ./generate_target_dnadiff_matrices.py ../fixtures/viral_example ../fixtures/viral_example/matrices
	cd tests/generate_fixtures/; ./generate_target_dnadiff_matrices.py ../fixtures/bad_alignments ../fixtures/bad_alignments/matrices
	cd tests/generate_fixtures/; ./generate_target_sourmash_files.py
	cd tests/generate_fixtures/; ./generate_target_sourmash_matrices.py
	@echo ""
	@echo "WARNING: If any tool version has changed, the generated test files"
	@echo "under fixtures/ may have changed (check using 'git diff'), and"
	@echo "thus the test output too. Beware!"
	@echo ""
	@echo "Please run generate_pyani_anib_matrices.sh manually with pyani v0.2"
	@echo "Please run generate_pyani_anim_matrices.sh manually with pyani v0.3"

# Run tests
# When the tests complete, the coverage output will be opened in a browser
# See also pyproject.toml setting addopts which adds --doctest-modules etc
test:
	@python -m pytest -n auto --cov-report=html --cov=pyani_plus -v && open htmlcov/index.html

# Clean up test output
clean_test:
	@rm -rf htmlcov
	@rm -rf tests/nucmer_filter_output
	@rm -rf tests/nucmer_delta_output
