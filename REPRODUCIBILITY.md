# Reproducibility

## Requirements

- Python 3.10 or newer
- `mpmath` 1.3.0 or a compatible later release
- `pdflatex` for rebuilding the PDF

Install the Python dependency with:

```bash
python -m pip install mpmath
```

## Run the numerical suite

From the repository root:

```bash
python verify_pi_g_identities.py --output-dir output
```

Expected overall verdict:

```text
PI_G_CROSS_DOMAIN_CONSISTENCY_CHECKS_PASS
```

The command writes:

```text
output/pi_g_verification_manifest.json
output/pi_g_verification_output.txt
output/pi_g_verification_results.csv
output/pi_g_verification_summary.txt
```

## Rebuild the report

```bash
python build_pi_g_report.py --compile
```

This command:

1. reruns the numerical suite;
2. synchronizes the construction and verification listings in the LaTeX source;
3. synchronizes the canonical output listing;
4. compiles the report twice with `pdflatex`.

To compile without rerunning the numerical checks:

```bash
python build_pi_g_report.py --skip-verification --compile
```

## Interpretation

The construction module is separated from the verification module. Classical
pi is not imported or referenced by `pi_g_construction.py`. The verification
module uses classical pi only after construction for a labelled diagnostic
comparison and for library-level consistency checks.

The suite is numerical corroboration, not a proof of equality.
