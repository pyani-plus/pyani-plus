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

import logging
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

from pyani_plus import PROGRESS_BAR_COLUMNS, db_orm, log_sys_exit, setup_logger
from pyani_plus.public_cli_args import ToolExecutor


class ShowProgress(str, Enum):
    """How to show progress of the workflow execution."""

    quiet = "quiet"
    bar = "bar"
    spin = "spin"


def progress_bar_via_db_comparisons(
    database: Path, run_id: int, interval: float = 1.0
) -> None:
    """Show a progress bar based on monitoring the DB entries.

    Will only self-terminate once all the run's comparisons are
    in the database! Up to the caller to terminate it early.
    """
    # Ignore any keyboard interrupte - let main thread handle it
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    logger = setup_logger(None)  # this is not on the main thead!
    session = db_orm.connect_to_db(logger, database)
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
    logger: logging.Logger,
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
    log: Path | None = None,
) -> None:
    """Run snakemake with a progress bar.

    The datatabase and run_id are only required with a progress bar,
    which will be used for live updates.

    In quiet or spinner mode the DB is only accessed except by the workflow
    itself, and need not be passed to this function.
    """
    msg = f"Preparing to call snakemake on '{workflow_name}'"
    logger.debug(msg)
    success = False
    if temp:
        params["temp"] = str(temp.resolve())
    if log:
        params["log"] = str(log.resolve())

    show_progress_bar = display == ShowProgress.bar
    if show_progress_bar and (database is None or run_id is None):
        msg = "Both database and run_id are required with display as progress bar"
        log_sys_exit(logger, msg)

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
        logger.debug("Calling snakemake without progress bar")
        success = args_to_api(args, parser)
    else:
        logger.debug("Calling snakemake with progress bar")
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
        log_sys_exit(logger, "Snakemake workflow failed")  # pragma: no cover
    logger.debug("Snakemake finished successfully")
