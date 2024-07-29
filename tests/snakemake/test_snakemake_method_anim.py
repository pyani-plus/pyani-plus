import shutil
import pytest

from pyani_plus.snakemake import snakemake_scheduler


@pytest.fixture
def config_filter_args(
    anim_nucmer_targets_filter_indir, anim_nucmer_targets_filter_values_outdir
):
    """Configuration settings for testing snakemake filter rule.

    We take the output directories for the MUMmer filter output and the
    small set of input genomes as arguments.
    """
    return {
        "outdir": str(anim_nucmer_targets_filter_values_outdir),
        "indir": str(anim_nucmer_targets_filter_indir),
        "cores": 8,
    }


def test_snakemake_rule_anim_method(
    anim_nucmer_targets_filter_values,
    config_filter_args,
):
    """Test nucmer filter snakemake wrapper

    Checks that the filter rule in the ANIm snakemake wrapper gives the
    expected output.

    If the output directory exists (i.e. the make clean_tests rule has not
    been run), the tests will automatically pass as snakemake will not
    attempt to re-run the rule. That would prevent us from seeing any
    introduced bugs, so we force re-running the rule by deleting the
    output directory before running the tests.
    """
    # # Remove the output directory to force re-running the snakemake rule
    # shutil.rmtree(anim_nucmer_targets_filter_outdir, ignore_errors=True)

    # Run snakemake wrapper
    runner = snakemake_scheduler.SnakemakeRunner("snakemake_anim_method.smk")

    runner.run_workflow(anim_nucmer_targets_filter_values, config_filter_args)
