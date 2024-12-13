The `viral_example` test set consists of three phage genomes that share similar regions, resulting in aligned regions. These viral genomes are currently used for testing `pyANI-plus` due to the short runtime of the test. This test set is used to check the reliability of `pyANI-plus`, ensuring correct parsing of intermediate files, accurate logging of values to the database, and the correct report of values in the output matrices.

This directory also includes the following additional subdirectories:
- `intermediates`: contains intermediate files used to test all ANI methods for the `viral_example` test case (e.g., `.filter` and `.delta` files for `ANIm`, `.sig` and `.csv`, for `sourmash` methods etc.).
- `matrices`: includes the expected `pyANI-plus` matrix outputs for the `viral_example` test case for each method implemented in `pyANI-plus`.
