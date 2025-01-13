# The MIT License
#
# Copyright (c) 2024-2025 University of Strathclyde
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
import signal
import sys
from enum import Enum
from pathlib import Path
from time import sleep

from rich.progress import Progress
from snakemake.cli import args_to_api, parse_args
from sqlalchemy import text
from sqlalchemy.exc import OperationalError

from pyani_plus import FASTA_EXTENSIONS, PROGRESS_BAR_COLUMNS, db_orm


class ToolExecutor(str, Enum):
    """Enum for the executor command line argument passed to snakemake."""

    local = "local"
    slurm = "slurm"


class ShowProgress(str, Enum):
    """How to show progress of the workflow execution."""

    quiet = "quiet"
    bar = "bar"
    spin = "spin"


def check_input_stems(indir: str) -> dict[str, Path]:
    """Check input files against approved list of extensions.

    If duplicate stems with approved extensions are present
    raise a ValueError.
    """
    extensions = tuple(FASTA_EXTENSIONS.union(_ + ".gz" for _ in FASTA_EXTENSIONS))

    stems = [_.stem for _ in Path(indir).glob("*") if _.name.endswith(extensions)]

    if len(stems) == len(set(stems)):
        input_files = {
            _.stem: _ for _ in Path(indir).glob("*") if _.name.endswith(extensions)
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
    # Ignore any keyboard interrupte - let main thread handle it
    signal.signal(signal.SIGINT, signal.SIG_IGN)

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
            try:
                db_version = (
                    session.connection().execute(text("PRAGMA data_version;")).one()[0]
                )
            except OperationalError:  # pragma: no cover
                # Probably another process is updating the DB, try again soon
                continue  # pragma: no cover
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
    display: ShowProgress = ShowProgress.quiet,
    database: Path | None = None,
    run_id: int | None = None,
    interval: float = 0.5,
    temp: Path | None = None,
) -> None:
    """Run snakemake with a progress bar.

    The datatabase and run_id are only required with a progress bar,
    which will be used for live updates.

    In quiet or spinner mode the DB is only accessed except by the workflow
    itself, and need not be passed to this function.
    """
    success = False
    # Including the --temp as part of the parameter as a way to avoid
    # adding --temp without an argument:
    if temp:
        params["temp"] = f"--temp {temp.resolve()}"
    else:
        # With empty string or '' or "" snakemake puts the literal string
        # None into the command line. So can't use that.
        # So as a horrible hack, have to put in something, so repeat --quiet
        params["temp"] = "--quiet"

    show_progress_bar = display == ShowProgress.bar
    if show_progress_bar and (database is None or run_id is None):
        msg = "Both database and run_id are required with display as progress bar"
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
            # No argument for snakemake quiet mode with progress bar or spinner,
            *(["rules"] if display == ShowProgress.quiet else []),
            "--executor",
            executor.value,
            "--directory",
            str(working_directory),
            "--snakefile",
            str(snakefile),
        ]
        + [str(_) for _ in targets]
    )
    args.config = [f"{k}='{v}'" for k, v in params.items()]
    if display == ShowProgress.quiet:
        success = args_to_api(args, parser)
    elif display == ShowProgress.spin:
        with Progress(*PROGRESS_BAR_COLUMNS) as progress:
            # Not quite the visual look I want, but close:
            task = progress.add_task("Comparing pairs", total=None)
            success = args_to_api(args, parser)
            progress.update(task, advance=len(targets), total=len(targets))
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
            daemon=True,
        )
        p.start()

        # Call snakemake! This seems to catch in KeyboardInterrupt itself
        success = args_to_api(args, parser)

        if p.is_alive():
            # Progress bar should have finished, perhaps final update pending...
            sleep(interval)
            p.terminate()
        p.join()

    if not success:
        # Writing a reliable test to trigger this has proved difficult,
        # so marking this as not expecting any test coverage.

        # Ensure exit message starts on a new line after interrupted progress bar
        print()  # pragma: no cover  # noqa: T201
        sys.exit("Snakemake workflow failed")  # pragma: no cover
