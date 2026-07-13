#!/usr/bin/env python3
"""Synchronize and optionally compile the Baker pi_g verification report."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
PAPER_DIR = ROOT / "paper"
TEX_PATH = PAPER_DIR / "pi_g_verification_report.tex"
OUTPUT_DIR = ROOT / "output"


def replace_listing(source: str, begin_line_pattern: str, body: str) -> str:
    pattern = re.compile(
        rf"({begin_line_pattern}\n).*?(\n\\end{{lstlisting}})",
        flags=re.DOTALL,
    )
    replacement = rf"\1{body.rstrip()}\2"
    updated, count = pattern.subn(replacement, source, count=1)
    if count != 1:
        raise RuntimeError(f"Could not uniquely replace listing: {begin_line_pattern}")
    return updated


def run_verification() -> None:
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "verify_pi_g_identities.py"),
            "--output-dir",
            str(OUTPUT_DIR),
        ],
        cwd=ROOT,
        check=True,
    )


def synchronize_report() -> None:
    source = TEX_PATH.read_text(encoding="utf-8")
    construction = (ROOT / "pi_g_construction.py").read_text(encoding="utf-8")
    verification = (ROOT / "verify_pi_g_identities.py").read_text(encoding="utf-8")
    output_log = (OUTPUT_DIR / "pi_g_verification_output.txt").read_text(
        encoding="utf-8"
    )

    source = replace_listing(
        source,
        r"\\begin\{lstlisting\}\[style=pythonstyle,caption=\{\\texttt\{pi\\_g\\_construction\.py\}\}\]",
        construction,
    )
    source = replace_listing(
        source,
        r"\\begin\{lstlisting\}\[style=pythonstyle,caption=\{\\texttt\{verify\\_pi\\_g\\_identities\.py\}\}\]",
        verification,
    )
    source = replace_listing(
        source,
        r"\\begin\{lstlisting\}\[style=outputstyle,caption=\{Canonical verification output\}\]",
        output_log,
    )
    TEX_PATH.write_text(source, encoding="utf-8")


def compile_report() -> None:
    command = [
        "pdflatex",
        "-interaction=nonstopmode",
        "-halt-on-error",
        TEX_PATH.name,
    ]
    for _ in range(2):
        subprocess.run(command, cwd=PAPER_DIR, check=True)

    for suffix in (".aux", ".log", ".out", ".toc"):
        path = TEX_PATH.with_suffix(suffix)
        if path.exists():
            path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compile",
        action="store_true",
        help="Compile the synchronized LaTeX source with pdflatex.",
    )
    parser.add_argument(
        "--skip-verification",
        action="store_true",
        help="Reuse the existing output ledger instead of rerunning the suite.",
    )
    args = parser.parse_args()

    if not args.skip_verification:
        run_verification()

    synchronize_report()

    if args.compile:
        compile_report()

    print(f"report_source: {TEX_PATH}")
    print(f"report_pdf: {TEX_PATH.with_suffix('.pdf')}")
    print(f"output_dir: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
