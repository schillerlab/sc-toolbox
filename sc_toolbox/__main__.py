#!/usr/bin/env python
"""Command-line interface."""
import logging
import os
import sys

import click
import rich.logging
from rich import print
from rich import traceback
from sc_toolbox.cli.commands.create import ProjectCreator
from sc_toolbox.cli.commands.upgrade import UpgradeCommand

WD = os.path.dirname(__file__)
log = logging.getLogger()


def main() -> None:
    traceback.install(width=200, word_wrap=True)
    print(
        r"""[bold blue]
  ██████  ▄████▄  ▄▄▄█████▓ ▒█████   ▒█████   ██▓     ▄▄▄▄    ▒█████  ▒██   ██▒
▒██    ▒ ▒██▀ ▀█  ▓  ██▒ ▓▒▒██▒  ██▒▒██▒  ██▒▓██▒    ▓█████▄ ▒██▒  ██▒▒▒ █ █ ▒░
░ ▓██▄   ▒▓█    ▄ ▒ ▓██░ ▒░▒██░  ██▒▒██░  ██▒▒██░    ▒██▒ ▄██▒██░  ██▒░░  █   ░
  ▒   ██▒▒▓▓▄ ▄██▒░ ▓██▓ ░ ▒██   ██░▒██   ██░▒██░    ▒██░█▀  ▒██   ██░ ░ █ █ ▒
▒██████▒▒▒ ▓███▀ ░  ▒██▒ ░ ░ ████▓▒░░ ████▓▒░░██████▒░▓█  ▀█▓░ ████▓▒░▒██▒ ▒██▒
▒ ▒▓▒ ▒ ░░ ░▒ ▒  ░  ▒ ░░   ░ ▒░▒░▒░ ░ ▒░▒░▒░ ░ ▒░▓  ░░▒▓███▀▒░ ▒░▒░▒░ ▒▒ ░ ░▓ ░
░ ░▒  ░ ░  ░  ▒       ░      ░ ▒ ▒░   ░ ▒ ▒░ ░ ░ ▒  ░▒░▒   ░   ░ ▒ ▒░ ░░   ░▒ ░
░  ░  ░  ░          ░      ░ ░ ░ ▒  ░ ░ ░ ▒    ░ ░    ░    ░ ░ ░ ░ ▒   ░    ░
      ░  ░ ░                   ░ ░      ░ ░      ░  ░ ░          ░ ░   ░    ░
         ░                                                 ░
    """
    )
    print("[bold blue]Run [green]sc-toolbox --help [blue]for an overview of all commands\n")

    # Is the latest sc-toolbox version installed? Upgrade if not!
    if not UpgradeCommand.check_sc_toolbox_latest():
        print("[bold blue]Run [green]sc-toolbox upgrade [blue]to get the latest version.")
    sc_toolbox_cli()  # type: ignore


@click.group()
@click.option("-v", "--verbose", is_flag=True, default=False, help="Enable verbose output (print debug statements).")
@click.option("-l", "--log-file", help="Save a verbose log to a file.")
def sc_toolbox_cli(verbose, log_file):
    """
    Create state of the art projects from production ready templates.
    """
    # Set the base logger to output DEBUG
    log.setLevel(logging.DEBUG)

    # Set up logs to the console
    log.addHandler(
        rich.logging.RichHandler(
            level=logging.DEBUG if verbose else logging.INFO,
            console=rich.console.Console(file=sys.stderr),
            show_time=True,
            markup=True,
        )
    )

    # Set up logs to a file if we asked for one
    if log_file:
        log_fh = logging.FileHandler(log_file, encoding="utf-8")
        log_fh.setLevel(logging.DEBUG)
        log_fh.setFormatter(logging.Formatter("[%(asctime)s] %(name)-20s [%(levelname)-7s]  %(message)s"))
        log.addHandler(log_fh)


@sc_toolbox_cli.command(short_help="Create a new project based on templates")
def create() -> None:
    """
    Create a new project based on an existing template.
    Usually includes a Docker container, a Conda environment and notebooks.
    """
    project_creator = ProjectCreator()
    project_creator.create_project()


@sc_toolbox_cli.command(short_help="Check for a newer version of sc-toolbox and upgrade if required.")  # type: ignore
def upgrade() -> None:
    """
    Checks whether the locally installed version of mlf-core is the latest.
    If not pip will be invoked to upgrade mlf-core to the latest version.
    """
    UpgradeCommand.check_upgrade_sc_toolbox()


if __name__ == "__main__":
    traceback.install()
    main(prog_name="sc-toolbox")  # type: ignore
