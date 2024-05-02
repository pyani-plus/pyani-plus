import shutil
import pytest

from pathlib import Path

from pyani_plus.snakemake import fastani


def compare_fastani_files(file1, file2):
    """Compare two fastANI files.

    This function expects two text files as input and returns True if the content
    of the files is the same, and False if the two files differ.

    As the Path to both Query and Reference might be diffrent,
    we will only consider file name.
    """

    with file1.open() as if1, file2.open() as if2:
        for line1, line2 in zip(if1.readlines(), if2.readlines(), strict=False):
            # Split lines into items
            line1 = line1.split("\t")
            line2 = line2.split("\t")

            # Extract file names from paths
            qry1 = line1[0].split("/")[-1]
            ref1 = line1[1].split("/")[-1]
            qry2 = line2[0].split("/")[-1]
            ref2 = line2[1].split("/")[-1]

            # Extract ANI value
            ANI1 = line1[2]
            ANI2 = line2[2]

            # Extract count of bidirectional fragment mappings
            fragMAP1 = line1[3]
            fragMAP2 = line2[3]

            # Extract total query framents
            qryFrag1 = line1[4]
            qryFrag2 = line2[4]

            # Reconstruct fastANI format
            fastani1 = "\t".join([qry1, ref1, ANI1, fragMAP1, qryFrag1])
            fastani2 = "\t".join([qry2, ref2, ANI2, fragMAP2, qryFrag2])
            if fastani1 != fastani2:
                return False
        return True


@pytest.fixture
def config_fastani_args(fastani_targets_outdir, input_genomes_small):
    """Configuration settings for testing snakemake fastANI rule."""
    return {
        "outdir": fastani_targets_outdir,
        "indir": input_genomes_small,
        "cores": 8,
        "fragLen": 3000,
        "kmerSize": 16,
        "minFrac": 0.2,
    }


def test_compare_fastani_files(
    fastani_targets,
    fastani_expected_targets,
    fastani_targets_outdir,
    config_fastani_args,
):
    """Test Snakemake fastani output.

    This pytest verifies that the fastani rule in the fastani Snakemake wrapper
    generates files with the expected content. The verification is
    achieved by comparing the output returned by the Snakemake workflow with
    the expected output.

    Before comparing the outputs, any existing content generated by the Snakemake
    filter rule is removed, and the rule is re-run to ensure a clean environment
    and to avoid any interference from previous executions.
    """

    # Remove the output directory to force re-running the snakemake rule
    shutil.rmtree(fastani_targets_outdir, ignore_errors=True)

    # Run snakemake wrapper
    fastani.run_workflow(fastani_targets, config_fastani_args)

    for fname in fastani_expected_targets.glob("*fastani"):
        for fname2 in fastani_targets_outdir.glob("*fastani"):
            if str(fname.stem) == str(fname2.stem):
                assert compare_fastani_files(fname, fname2) is True
