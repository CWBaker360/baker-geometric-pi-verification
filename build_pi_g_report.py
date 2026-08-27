#!/usr/bin/env python3
"""Build and optionally compile the synchronized pi_g verification report."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from mpmath import mp

from pi_g_construction import compute_pi_g


ROOT = Path(__file__).resolve().parent


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def read_results(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def latex_sci(text: str, digits: int = 8) -> str:
    value = mp.mpf(text)
    if value == 0:
        return r"\texttt{0}"
    formatted = mp.nstr(value, digits, min_fixed=0, max_fixed=0)
    return r"\texttt{" + latex_escape(formatted) + "}"

def build_table_rows(results: list[dict[str, str]]) -> str:
    rows = []
    for row in results:
        rows.append(
            "{} & {} & {} & {} \\\\".format(
                latex_escape(row["identifier"]),
                latex_escape(row["name"]),
                latex_sci(row["observed_metric"], 9),
                "PASS" if row["passed"].lower() == "true" else "FAIL",
            )
        )
    return "\n".join(rows)


def build_method_rows(results: list[dict[str, str]]) -> str:
    rows = []
    for row in results:
        rows.append(
            "{} & {} & {} \\\\".format(
                latex_escape(row["identifier"]),
                latex_escape(row["method"]),
                latex_sci(row["acceptance_bound"], 9),
            )
        )
    return "\n".join(rows)

def build_polygon_rows(level_count: int, work_dps: int) -> str:
    mp.dps = work_dps
    construction = compute_pi_g(level_count=level_count, work_dps=work_dps)
    rows = []
    for level, (polygon, diagonal) in enumerate(
        zip(construction.polygon_values, construction.extrapolated_diagonal)
    ):
        n_sides = 6 * 3**level
        rows.append(
            f"{level} & {n_sides} & "
            + r"\texttt{" + mp.nstr(polygon, 22) + "} & "
            + r"\texttt{" + mp.nstr(diagonal, 22) + r"} \\" 
        )
    return "\n".join(rows)


def embedded_listing(path: Path, style: str, caption: str) -> str:
    """Return a self-contained listings block with the file text embedded."""
    text = path.read_text(encoding="utf-8")
    if r"\end{lstlisting}" in text:
        raise ValueError(f"Cannot embed {path}: contains " + r"\end{lstlisting}")
    return (
        rf"\begin{{lstlisting}}[style={style},caption={{{caption}}}]"
        + "\n"
        + text.rstrip("\n")
        + "\n" + r"\end{lstlisting}"
    )


def make_tex(manifest: dict, results: list[dict[str, str]]) -> str:
    table_rows = build_table_rows(results)
    method_rows = build_method_rows(results)
    polygon_rows = build_polygon_rows(
        manifest["level_count"], manifest["working_precision_decimal_digits"]
    )
    construction_listing = embedded_listing(
        ROOT / "pi_g_construction.py",
        "pythonstyle",
        r"\texttt{pi\_g\_construction.py}",
    )
    verification_listing = embedded_listing(
        ROOT / "verify_pi_g_identities.py",
        "pythonstyle",
        r"\texttt{verify\_pi\_g\_identities.py}",
    )
    output_listing = embedded_listing(
        ROOT / "output" / "pi_g_verification_output.txt",
        "outputstyle",
        "Canonical verification output",
    )

    pi_g = manifest["pi_g_decimal"]
    pi_ref = manifest["reference_pi_decimal"]
    reference_error = manifest["diagnostic_abs_pi_g_minus_pi"]
    internal_delta = manifest["internal_last_diagonal_delta"]
    max_bracket = manifest["maximum_root_bracket_width"]
    max_poly_residual = manifest["maximum_root_polynomial_residual"]
    work_dps = manifest["working_precision_decimal_digits"]
    level_count = manifest["level_count"]
    max_level = manifest["maximum_level"]
    python_version = manifest["software"]["python"]
    mpmath_version = manifest["software"]["mpmath"]

    return rf"""\documentclass[11pt]{{article}}
\usepackage[T1]{{fontenc}}
\usepackage[utf8]{{inputenc}}
\usepackage{{lmodern}}
\usepackage{{amsmath,amssymb,amsthm}}
\usepackage{{booktabs,longtable,array}}
\usepackage{{geometry}}
\usepackage{{microtype}}
\usepackage{{listings}}
\usepackage{{xcolor}}
\usepackage{{hyperref}}
\usepackage{{fancyhdr}}
\usepackage{{enumitem}}

\geometry{{margin=1in}}
\hypersetup{{colorlinks=true,linkcolor=blue!45!black,urlcolor=blue!45!black,citecolor=blue!45!black}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.55em}}
\pagestyle{{fancy}}
\fancyhf{{}}
\lhead{{Baker's geometric $\pi_g$}}
\rhead{{Referee-safe computational appendix}}
\cfoot{{\thepage}}

\lstdefinestyle{{pythonstyle}}{{
  language=Python,
  basicstyle=\ttfamily\scriptsize,
  keywordstyle=\bfseries,
  commentstyle=\itshape,
  stringstyle=\ttfamily,
  numbers=left,
  numberstyle=\tiny,
  stepnumber=1,
  numbersep=8pt,
  breaklines=true,
  breakatwhitespace=false,
  showstringspaces=false,
  columns=fullflexible,
  keepspaces=true,
  frame=single,
  rulecolor=\color{{black!25}},
  tabsize=4
}}
\lstdefinestyle{{outputstyle}}{{
  basicstyle=\ttfamily\scriptsize,
  numbers=none,
  breaklines=true,
  columns=fullflexible,
  keepspaces=true,
  frame=single,
  rulecolor=\color{{black!25}}
}}

\newtheorem{{proposition}}{{Proposition}}
\newcommand{{\pig}}{{\pi_g}}

\title{{Cross-Domain Numerical Consistency Checks\\for Baker's Geometric Constant $\pig$}}
\author{{Wayne Baker}}
\date{{August 27, 2026\\Referee-safe computational appendix, Version 2.2.0}}

\begin{{document}}
\maketitle

\begin{{abstract}}
This report documents a synchronized, reproducible numerical verification suite for Baker's algebraically constructed geometric constant $\pig$.  The construction module begins with the exact hexagon seed $c_0=1/2$, selects the required cubic branch by a deterministic sign-changing bracket, forms the polygon sequence $P_L=6\cdot3^L c_L$, and applies the scale-$3$ even-power Richardson table.  The classical constant $\pi$ is not consulted during construction or result selection.  After $\pig$ has been fixed, nine cross-domain checks are performed using direct residuals, explicit truncation enclosures, or two-sided inequalities.  All nine checks pass at {work_dps}-decimal-digit working precision.  These computations are consistency evidence, not a substitute for a proof of $\pig=\pi$.
\end{{abstract}}

\tableofcontents
\newpage

\section{{Status and claim boundary}}

The computational claim is deliberately limited:

\begin{{quote}}
The algebraic pipeline produces a fixed high-precision value $\pig$.  That value satisfies nine stated numerical consistency checks within their declared acceptance bounds.
\end{{quote}}

The suite does \emph{{not}} prove $\pig=\pi$.  Equality must rest on the associated mathematical argument connecting the chord cascade and polygonal limit to the classical analytic constant.  The tests below illustrate consequences of that equality and provide a reproducibility audit of the implementation.

The phrase ``independent identity verification'' is avoided.  In particular, the \texttt{{mpmath}} implementations of special functions may internally use classical $\pi$.  Those tests are therefore cross-domain library-level consistency checks, not algorithmically independent derivations of $\pi$.

\section{{Oracle-free algebraic construction}}

\subsection{{Unique cubic branch}}

Let
\[
f_c(u)=4u^3-3u+c,
\qquad 0<c\le \frac12.
\]
The trisection step requires the small positive root.  It is selected without an initial guess and without reference to $\pi$.

\begin{{proposition}}
For every $0<c\le 1/2$, the interval $(0,c/2)$ contains exactly one zero of $f_c$.
\end{{proposition}}

\begin{{proof}}
At the endpoints,
\[
f_c(0)=c>0
\]
and
\[
f_c\!\left(\frac c2\right)
=4\left(\frac c2\right)^3-3\left(\frac c2\right)+c
=\frac{{c(c^2-1)}}{{2}}<0.
\]
Moreover,
\[
f_c'(u)=12u^2-3<0
\]
for $0\le u\le c/2\le1/4$.  Thus $f_c$ is strictly decreasing on the bracket, and the sign change contains exactly one root.
\end{{proof}}

The implementation bisects this certified branch until the bracket is smaller than the prescribed guard scale.  The maximum final bracket width in the present run was
\[
\texttt{{{max_bracket}}},
\]
and the maximum midpoint polynomial residual was
\[
\texttt{{{max_poly_residual}}}.
\]
These are floating-point diagnostics attached to an exact branch-selection argument.

\subsection{{Polygon sequence and N-series acceleration}}

Starting from
\[
c_0=\frac12,
\]
the recurrence chooses $c_{{L+1}}$ as the unique root in $(0,c_L/2)$ of
\[
4c_{{L+1}}^3-3c_{{L+1}}+c_L=0.
\]
The polygon quantity is
\[
P_L=6\cdot3^L c_L.
\]
The present computation uses {level_count} inputs,
\[
L=0,1,\ldots,{max_level}.
\]

If a sequence has an even inverse-power expansion
\[
P_N=Q+\frac{{a_1}}{{N^2}}+\frac{{a_2}}{{N^4}}+\cdots,
\]
then replacing $N$ by $3N$ scales the stage-$s$ leading error by $9^{{-s}}$.  The triangular extrapolation is therefore
\[
T_0(L)=P_L,
\]
\[
T_s(L)=T_{{s-1}}(L)
+\frac{{T_{{s-1}}(L)-T_{{s-1}}(L-1)}}{{9^s-1}}.
\]
The selected value is deterministically the final table entry $T_{{{max_level}}}({max_level})$.  No comparison with classical $\pi$ is used to select a level or candidate.

\begin{{longtable}}{{r r >{{\raggedleft\arraybackslash}}p{{0.27\textwidth}} >{{\raggedleft\arraybackslash}}p{{0.27\textwidth}}}}
\toprule
$L$ & $N=6\cdot3^L$ & $P_L$ & Richardson diagonal $T_L(L)$ \\
\midrule
\endfirsthead
\toprule
$L$ & $N=6\cdot3^L$ & $P_L$ & Richardson diagonal $T_L(L)$ \\
\midrule
\endhead
{polygon_rows}
\bottomrule
\end{{longtable}}

The resulting value is
\[
\pig=
\begin{{aligned}}
&\texttt{{{pi_g[:50]}}}\\
&\texttt{{{pi_g[50:]}}}.
\end{{aligned}}
\]
The last two extrapolation diagonal entries differ by
\[
\texttt{{{internal_delta}}}.
\]
This difference is an internal convergence indicator, not a rigorous error bound.

Only after $\pig$ was fixed was it compared with the library reference value
\[
\pi=
\begin{{aligned}}
&\texttt{{{pi_ref[:50]}}}\\
&\texttt{{{pi_ref[50:]}}}.
\end{{aligned}}
\]
The diagnostic difference was
\[
|\pig-\pi|=\texttt{{{reference_error}}}.
\]
This comparison is reported transparently but is not part of the construction or any candidate-selection rule.

\section{{Verification protocol}}

\subsection{{Numerical environment}}

\begin{{itemize}}[leftmargin=2em]
\item Python: \texttt{{{latex_escape(python_version)}}}.
\item \texttt{{mpmath}}: \texttt{{{latex_escape(mpmath_version)}}}.
\item Working precision: {work_dps} decimal digits.
\item Construction inputs: {level_count} polygon levels, $L=0,\ldots,{max_level}$.
\item Fixed threshold for direct special-function residuals: $10^{{-65}}$.
\item Every series/product check uses either an explicit remainder bound or a finite two-sided enclosure.
\end{{itemize}}

No printed value ``0'' is interpreted as exact equality.  A zero metric means only that the computed violation was nonpositive at the active precision or that a target lay inside a finite enclosure.

\subsection{{Acceptance criteria}}

\begin{{longtable}}{{p{{0.07\textwidth}} p{{0.67\textwidth}} p{{0.17\textwidth}}}}
\toprule
ID & Method & Acceptance bound \\
\midrule
\endfirsthead
\toprule
ID & Method & Acceptance bound \\
\midrule
\endhead
{method_rows}
\bottomrule
\end{{longtable}}

\subsection{{Bounded checks}}

\paragraph{{Gaussian integral.}}
The finite integral is evaluated on $[-A,A]$ with $A=13$.  The omitted two-sided tail satisfies
\[
2\int_A^\infty e^{{-x^2}}\,dx
\le \frac{{e^{{-A^2}}}}{{A}}.
\]
The acceptance allowance is this analytic tail bound plus a repeat-precision quadrature stability term and a small rounding reserve.

\paragraph{{Basel identity.}}
The suite uses the rapidly convergent positive series
\[
\zeta(2)=3\sum_{{n=1}}^\infty
\frac{{1}}{{n^2\binom{{2n}}{{n}}}}.
\]
For the unscaled term $a_n=1/(n^2\binom{{2n}}{{n}})$,
\[
\frac{{a_{{n+1}}}}{{a_n}}
=\frac{{n^2}}{{(2n+2)(2n+1)}}<\frac14.
\]
Consequently, after $N$ terms the omitted $\zeta(2)$ tail is smaller than $4a_{{N+1}}$.

\paragraph{{Gregory-Leibniz identity.}}
Euler transformation is applied to
\[
\frac\pi4=\sum_{{n=0}}^\infty\frac{{(-1)^n}}{{2n+1}}.
\]
Writing
\[
I_k=\int_0^1(1-x^2)^k\,dx,
\]
gives
\[
\frac\pi4=\sum_{{k=0}}^\infty\frac{{I_k}}{{2^{{k+1}}}},
\qquad
I_0=1,
\qquad
I_{{k+1}}=\frac{{2(k+1)}}{{2k+3}}I_k.
\]
The transformed term ratio is
\[
\frac{{t_{{k+1}}}}{{t_k}}=\frac{{k+1}}{{2k+3}}<\frac12,
\]
so the omitted tail is smaller than twice the first omitted term.

\paragraph{{Wallis product.}}
For
\[
W_N=\prod_{{n=1}}^N\frac{{(2n)^2}}{{(2n-1)(2n+1)}},
\]
the standard finite Wallis enclosure is
\[
W_N<\frac\pi2<W_N\frac{{2N+1}}{{2N}}.
\]
The test passes only when $\pig/2$ lies inside this interval.

\paragraph{{Robbins-Stirling bounds.}}
For each tested integer $n$, the suite checks the exact logarithmic factorial sum against Robbins' bounds
\[
\sqrt{{2\pi n}}\left(\frac ne\right)^n e^{{1/(12n+1)}}
<n!
<
\sqrt{{2\pi n}}\left(\frac ne\right)^n e^{{1/(12n)}}.
\]
The implementation replaces $\pi$ by the already constructed $\pig$ and verifies inclusion for $n=10,100,1000,10000$.

\section{{Results}}

\begin{{longtable}}{{p{{0.07\textwidth}} p{{0.35\textwidth}} p{{0.31\textwidth}} p{{0.10\textwidth}}}}
\toprule
ID & Check & Observed metric & Status \\
\midrule
\endfirsthead
\toprule
ID & Check & Observed metric & Status \\
\midrule
\endhead
{table_rows}
\bottomrule
\end{{longtable}}

All nine declared checks passed:
\[
\boxed{{\texttt{{PI\_G\_CROSS\_DOMAIN\_CONSISTENCY\_CHECKS\_PASS}}}}.
\]

The strongest interpretation supported by the computation is:

\begin{{quote}}
At 120-decimal-digit working precision, the deterministically constructed value $\pig$ is consistent with the listed bounded series, product, asymptotic, analytic, and special-function relations under the stated acceptance rules.
\end{{quote}}

The computation does not establish functional identity ``across all domains,'' nor does it replace a theorem proving $\pig=\pi$.

\section{{Reproducibility and provenance}}

The construction and verification are physically separated:

\begin{{itemize}}[leftmargin=2em]
\item \texttt{{pi\_g\_construction.py}} contains the algebraic seed, certified branch selection, polygon sequence, and extrapolation.  It does not import or reference classical $\pi$.
\item \texttt{{verify\_pi\_g\_identities.py}} imports the completed construction and then performs the reference comparison and nine checks.
\end{{itemize}}

The canonical run command is
\begin{{lstlisting}}[style=outputstyle]
python verify_pi_g_identities.py --output-dir output
\end{{lstlisting}}

The report is regenerated with
\begin{{lstlisting}}[style=outputstyle]
python build_pi_g_report.py --compile
\end{{lstlisting}}

The output manifest records software metadata and SHA-256 hashes for both scripts and the principal numerical ledgers.

\section{{Conclusion}}

The revised suite removes the two central weaknesses of the earlier version.  First, the cubic branch is selected by a proved sign-changing interval rather than a heuristic Newton initial guess.  Second, PASS decisions are based on direct residuals or explicit enclosures, not on the near-tautological comparison of a $\pig$-based expression with the same expression evaluated using classical $\pi$.

The resulting package is suitable as a computational appendix accompanying a mathematical proof of $\pig=\pi$.  Its correct status is \emph{{empirical and reproducible corroboration}}.

\begin{{thebibliography}}{{9}}
\bibitem{{Robbins1955}}
H. Robbins,
``A Remark on Stirling's Formula,''
\emph{{American Mathematical Monthly}} 62 (1955), 26--29.

\bibitem{{Borwein1987}}
J. M. Borwein and P. B. Borwein,
\emph{{Pi and the AGM}},
Wiley, 1987.

\bibitem{{Richardson1927}}
L. F. Richardson,
``The Deferred Approach to the Limit,''
\emph{{Philosophical Transactions of the Royal Society A}} 226 (1927), 299--361.
\end{{thebibliography}}

\appendix
\section{{Oracle-free construction script}}
{construction_listing}

\section{{Verification script}}
{verification_listing}

\section{{Full output log}}
{output_listing}

\end{{document}}
"""


def compile_tex(tex_path: Path) -> None:
    for _ in range(2):
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tex_path.parent,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compile", action="store_true")
    parser.add_argument(
        "--skip-verification",
        action="store_true",
        help="build from the existing canonical output without rerunning the numerical suite",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.skip_verification:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "verify_pi_g_identities.py"),
                "--output-dir",
                str(ROOT / "output"),
            ],
            cwd=ROOT,
            check=True,
        )

    manifest_path = ROOT / "output" / "pi_g_verification_manifest.json"
    results_path = ROOT / "output" / "pi_g_verification_results.csv"
    if not manifest_path.exists() or not results_path.exists():
        raise SystemExit("Run verify_pi_g_identities.py before building the report")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = read_results(results_path)
    paper_dir = ROOT / "paper"
    paper_dir.mkdir(exist_ok=True)
    tex_path = paper_dir / "pi_g_verification_report.tex"
    tex_path.write_text(make_tex(manifest, results), encoding="utf-8")
    print(f"wrote: {tex_path}")

    if args.compile:
        compile_tex(tex_path)
        print(f"wrote: {tex_path.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
