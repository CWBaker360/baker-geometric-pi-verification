# Cross-Domain Numerical Consistency Checks for Baker's Geometric Constant pi_g

**Author:** Wayne Baker  
**Date:** August 27, 2026  
**Version:** 2.2.0  
**Status:** Referee-safe computational appendix / reproducibility archive

This repository contains the paper, algebraic construction module, verification
suite, numerical ledgers, and report builder for:

> **Cross-Domain Numerical Consistency Checks for Baker's Geometric Constant pi_g**

## Read the report

- [`paper/pi_g_verification_report.pdf`](paper/pi_g_verification_report.pdf)
- [`paper/pi_g_verification_report.tex`](paper/pi_g_verification_report.tex)

## Construction

The construction begins with the exact hexagon half-chord seed

\[
c_0=\frac12,
\]

then selects the unique small positive root of

\[
4u^3-3u+c=0
\]

inside the certified interval \((0,c/2)\). The polygon sequence is

\[
P_L=6\cdot 3^L c_L.
\]

A scale-3 even-power Richardson table is then applied, with the final entry
selected deterministically. Classical pi is not consulted during construction
or candidate selection.

## Verification suite

After pi_g is fixed, nine declared checks are evaluated:

- Gaussian integral with an analytic tail enclosure
- central-binomial Basel series with an explicit tail bound
- Euler-transformed Gregory-Leibniz series
- finite Wallis enclosure
- Robbins-Stirling two-sided bounds
- Gamma at one half
- Euler's identity
- Euler's reflection formula at three rational arguments
- Richardson-tableau a priori enclosure using Taylor remainders, fixed data-independent weights, and propagated cubic-root bracket width

All nine checks pass at 120-decimal-digit working precision.

## Claim boundary

These computations are reproducible numerical consistency evidence. They do
not prove that pi_g equals pi. The special-function checks are library-level
cross-domain checks and are not asserted to be algorithmically independent of
the numerical library's implementation of pi.

## Run the suite

From the repository root:

```bash
python verify_pi_g_identities.py --output-dir output
```

## Rebuild the report

```bash
python build_pi_g_report.py --compile
```

The builder reruns the verification suite, synchronizes the two Python listings
and canonical output log embedded in the LaTeX report, and compiles the PDF.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── REPRODUCIBILITY.md
├── LICENSE_CODE
├── LICENSE_CONTENT.md
├── LICENSE_NOTICE.md
├── CHANGELOG.md
├── REQUIREMENTS.txt
├── PACKAGE_MANIFEST.json
├── SHA256SUMS.txt
├── .gitignore
├── pi_g_construction.py
├── verify_pi_g_identities.py
├── build_pi_g_report.py
├── paper/
│   ├── pi_g_verification_report.tex
│   ├── pi_g_verification_report.pdf
│   └── README.md
├── output/
│   ├── pi_g_verification_manifest.json
│   ├── pi_g_verification_output.txt
│   ├── pi_g_verification_results.csv
│   ├── pi_g_verification_summary.txt
│   └── README.md
├── scripts/
│   └── README.md
└── docs/
    ├── abstract.md
    ├── claim_boundary.md
    ├── repository_description.md
    └── github_upload_checklist.md
```

## License

Python code is licensed under the MIT License. The paper, LaTeX source, documentation, and other non-code content are licensed under CC BY 4.0. See `LICENSE_NOTICE.md` for details.
