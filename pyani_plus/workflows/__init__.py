# The MIT License
#
# Copyright (c) 2024 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Snakemake workflows for ANI pairwise comparisons."""

import multiprocessing
import sys
from enum import Enum
from pathlib import Path
from time import sleep

from rich.progress import Progress
from snakemake.cli import args_to_api, parse_args
from sqlalchemy import text

from pyani_plus import FASTA_EXTENSIONS, PROGRESS_BAR_COLUMNS, db_orm


class ToolExecutor(str, Enum):
    """Enum for the executor command line argument passed to snakemake."""

    local = "local"
    slurm = "slurm"


def check_input_stems(indir: str) -> dict[str, Path]:
    """Check input files against approved list of extensions.

    If duplicate stems with approved extensions are present
    raise a ValueError.
    """
    stems = [_.stem for _ in Path(indir).glob("*") if _.suffix in FASTA_EXTENSIONS]

    if len(stems) == len(set(stems)):
        input_files = {
            _.stem: _ for _ in Path(indir).glob("*") if _.suffix in FASTA_EXTENSIONS
        }
    else:
        duplicates = [
            item for item in stems if stems.count(item) > 1 and item in set(stems)
        ]
        msg = (
            f"Duplicated stems found for {sorted(set(duplicates))}. Please investigate."
        )
        raise ValueError(msg)

    return input_files


def progress_bar_via_db_comparisons(
    database: Path, run_id: int, interval: float = 1.0
) -> None:
    """Show a progress bar based on monitoring the DB entries.

    Will only self-terminate once all the run's comparisons are
    in the database! Up to the caller to terminate it early.
    """
    session = db_orm.connect_to_db(database)
    run = session.query(db_orm.Run).where(db_orm.Run.run_id == run_id).one()

    already_done = run.comparisons().count()
    total = run.genomes.count() ** 2 - already_done

    done = 0
    old_db_version = -1
    with Progress(*PROGRESS_BAR_COLUMNS) as progress:
        task = progress.add_task("Comparing pairs", total=total)
        while done < total:
            sleep(interval)
            # Have there been any DB changes?
            db_version = (
                session.connection().execute(text("PRAGMA data_version;")).one()[0]
            )
            if old_db_version == db_version:
                continue
            old_db_version = db_version
            new = run.comparisons().count() - already_done - done
            progress.update(task, advance=new)
            done += new
    session.close()


def run_snakemake_with_progress_bar(  # noqa: PLR0913
    executor: ToolExecutor,
    workflow_name: str,
    targets: list[Path] | list[str],
    params: dict,
    working_directory: Path,
    *,
    show_progress_bar: bool = False,
    database: Path | None = None,
    run_id: int | None = None,
    interval: float = 0.5,
) -> None:
    """Run snakemake with a progress bar.

    The datatabase and run_id are only required with a progress bar.
    """
    if show_progress_bar and (database is None or run_id is None):
        msg = "Both database and run_id are required with show_progress_bar=True"
        raise ValueError(msg)

    # Check this now to give clear up-front error - this function will be
    # call from within the actual snakemake workflow files themselves:
    if "indir" in params:
        check_input_stems(params["indir"])

    # Path to anim snakemake file
    snakefile = Path(__file__).with_name(workflow_name)

    # Want snakemake to pull everything else from its profiles,
    # most easily controlled via setting $SNAKEMAKE_PROFILE
    parser, args = parse_args(
        [
            "--quiet",
            "all" if show_progress_bar else "rules",
            "--executor",
            executor.value,
            "--directory",
            str(working_directory),
            "--snakefile",
            str(snakefile),
        ]
        + [str(_) for _ in targets]
    )
    args.config = [f"{k}={v}" for k, v in params.items()]
    if not show_progress_bar:
        success = args_to_api(args, parser)
    else:
        # As of Python 3.8 onwards, the default on macOS ("Darwin") is "spawn"
        # As of Python 3.12, the default of "fork" on Linux triggers a deprecation warning.
        # This should match the defaults on Python 3.14 onwards.
        # Note mypy currently can't handle this dynamic situation, their issue #8603
        p = multiprocessing.get_context(  # type:ignore [attr-defined]
            "spawn" if sys.platform == "darwin" else "forkserver"
        ).Process(
            target=progress_bar_via_db_comparisons,
            args=(
                database,
                run_id,
                interval,
            ),
        )
        p.start()

        # Call snakemake!
        try:
            success = args_to_api(args, parser)
        except Exception as err:  # noqa: BLE001
            p.terminate()
            # Ensure traceback starts on a new line after interrupted progress bar
            print()  # noqa: T201
            raise err from None

        if p.is_alive():
            # Should have finished, but perhaps final update is pending...
            sleep(interval)
            p.terminate()
        p.join()

    if not success:
        sys.exit("Snakemake workflow failed")
