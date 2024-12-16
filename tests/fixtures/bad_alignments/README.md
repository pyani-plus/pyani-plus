The `bad_alignments` test set consists of two highly divergent phage genomes (`MGV-GENOME-0264574.fas` and `MGV-GENOME-0357962.fna`) with no shared regions, resulting in no alignments. This test set ensures that `pyANI-plus` correctly detects and handles such comparisons without errors, recording them as `NULL` in the database and output matrices.

This directory also includes the following additional subdirectories:
- `intermediates`: contains intermediate files used to test all ANI methods for the `bad_alignments` test case (e.g., `.filter` and `.delta` files for `ANIm`, `.sig` and `.csv`, for `sourmash` methods etc.).
- `matrices`: includes the expected `pyANI-plus` matrix outputs for the `bad_alignments` test case for each method implemented in `pyANI-plus`.
