# Changelog

## Version 2.2.0 — 2026-08-27

- Added check N9: a proved a priori enclosure for the final N-series Richardson diagonal.
- N9 combines the sine Taylor alternating-series remainder bound, the tableau's fixed data-independent linear weights, and propagated cubic-root bracket widths.
- The internal last-diagonal delta is retained as a diagnostic only, not a rigorous error bound.
- Verification count updated from eight to nine; canonical result is 9/9 PASS.
- Regenerated the synchronized report and repository metadata for Version 2.2.


## Version 2.1 — 2026-07-12

- Replaced heuristic cubic-root selection with a proved sign-changing bracket.
- Separated the oracle-free construction module from the verification module.
- Replaced near-tautological comparisons with direct residuals, explicit
  truncation bounds, or finite two-sided enclosures.
- Recorded eight passing cross-domain consistency checks at 120 decimal digits.
- Added the full construction source, verification source, and canonical output
  log as appendices.
- Added a synchronized report builder and repository-level reproducibility files.
