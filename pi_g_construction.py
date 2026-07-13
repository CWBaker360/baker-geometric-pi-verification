#!/usr/bin/env python3
"""Oracle-free construction kernel for Baker's geometric constant pi_g.

This module deliberately does not import or reference the classical constant pi.
It starts from the algebraic hexagon seed c_0 = 1/2, repeatedly selects the
unique small positive root of

    4 u^3 - 3 u + c = 0,

and applies the N-series Richardson table with scale factor 3.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from mpmath import mp


@dataclass(frozen=True)
class RootCertificate:
    """Numerical record for one deterministically bracketed cubic root."""

    value: mp.mpf
    lower: mp.mpf
    upper: mp.mpf
    width: mp.mpf
    iterations: int
    polynomial_residual: mp.mpf


@dataclass(frozen=True)
class PiGConstructionResult:
    """Complete output of the algebraic polygon/extrapolation pipeline."""

    value: mp.mpf
    level_count: int
    maximum_level: int
    work_dps: int
    polygon_values: tuple[mp.mpf, ...]
    extrapolated_diagonal: tuple[mp.mpf, ...]
    internal_last_diagonal_delta: mp.mpf
    maximum_root_bracket_width: mp.mpf
    maximum_root_polynomial_residual: mp.mpf
    root_iterations: tuple[int, ...]


def _cubic(u: mp.mpf, c: mp.mpf) -> mp.mpf:
    return 4 * u**3 - 3 * u + c


def small_positive_trisection_root(c: mp.mpf) -> RootCertificate:
    """Return the unique root in (0, c/2) by deterministic bisection.

    For 0 < c <= 1/2,

        f_c(0) = c > 0,
        f_c(c/2) = c(c^2 - 1)/2 < 0,

    and f'_c(u) = 12u^2 - 3 < 0 on [0, c/2].  Hence this interval
    contains exactly one root, which is the required small positive branch.

    The mathematical branch selection is exact.  The returned decimal value
    is a high-precision midpoint of the final sign-changing bracket.
    """

    c = mp.mpf(c)
    if not (mp.mpf("0") < c <= mp.mpf("0.5")):
        raise ValueError(f"Expected 0 < c <= 1/2, received {c!r}")

    lower = mp.mpf("0")
    upper = c / 2
    f_lower = _cubic(lower, c)
    f_upper = _cubic(upper, c)

    if not (f_lower > 0 and f_upper < 0):
        raise ArithmeticError("Initial cubic bracket does not change sign")

    # Preserve roughly 20 decimal guard digits inside the active mp context.
    target_width = mp.power(10, -(mp.dps - 20))
    maximum_iterations = max(256, 8 * mp.dps)

    iterations = 0
    for iterations in range(1, maximum_iterations + 1):
        midpoint = (lower + upper) / 2
        if midpoint == lower or midpoint == upper:
            break

        f_midpoint = _cubic(midpoint, c)
        if f_midpoint > 0:
            lower = midpoint
        elif f_midpoint < 0:
            upper = midpoint
        else:
            lower = upper = midpoint
            break

        if upper - lower <= target_width:
            break
    else:
        raise ArithmeticError("Bisection did not reach the requested bracket width")

    value = (lower + upper) / 2
    width = upper - lower
    residual = abs(_cubic(value, c))

    if width > target_width * 2:
        raise ArithmeticError(
            f"Final root bracket is wider than expected: {mp.nstr(width, 12)}"
        )

    return RootCertificate(
        value=value,
        lower=lower,
        upper=upper,
        width=width,
        iterations=iterations,
        polynomial_residual=residual,
    )


def richardson_n_series(values: Sequence[mp.mpf]) -> tuple[mp.mpf, tuple[mp.mpf, ...]]:
    """Apply the scale-3 even-power Richardson/N-series table.

    If the polygon error has an expansion in N^{-2}, N^{-4}, ..., replacing
    N by 3N multiplies the leading term at stage s by 9^{-s}.  Therefore the
    stage-s denominator is 9^s - 1.
    """

    if not values:
        raise ValueError("At least one polygon value is required")

    n = len(values)
    table: list[list[mp.mpf | None]] = [[mp.mpf(v) for v in values]]

    for stage in range(1, n):
        previous = table[-1]
        current: list[mp.mpf | None] = [None] * n
        denominator = mp.power(9, stage) - 1
        for column in range(stage, n):
            right = previous[column]
            left = previous[column - 1]
            if right is None or left is None:
                raise ArithmeticError("Malformed Richardson table")
            current[column] = right + (right - left) / denominator
        table.append(current)

    diagonal: list[mp.mpf] = []
    for stage in range(n):
        value = table[stage][stage]
        if value is None:
            raise ArithmeticError("Missing Richardson diagonal entry")
        diagonal.append(value)

    final = table[-1][-1]
    if final is None:
        raise ArithmeticError("Missing final Richardson value")

    return final, tuple(diagonal)


def compute_pi_g(*, level_count: int = 11, work_dps: int = 120) -> PiGConstructionResult:
    """Compute pi_g without consulting the classical constant pi.

    ``level_count=11`` means the polygon inputs L=0,...,10 are used.  The
    construction result is selected deterministically as the final entry of
    the corresponding Richardson table; no reference constant is consulted.
    """

    if level_count < 2:
        raise ValueError("level_count must be at least 2")
    if work_dps < 50:
        raise ValueError("work_dps must be at least 50")

    with mp.workdps(work_dps):
        c = mp.mpf(1) / 2
        polygon_values: list[mp.mpf] = []
        root_certificates: list[RootCertificate] = []

        for level in range(level_count):
            n_sides = mp.mpf(6) * mp.power(3, level)
            polygon_values.append(n_sides * c)

            if level + 1 < level_count:
                certificate = small_positive_trisection_root(c)
                root_certificates.append(certificate)
                c = certificate.value

        final, diagonal = richardson_n_series(polygon_values)
        internal_delta = abs(diagonal[-1] - diagonal[-2])

        maximum_width = max(cert.width for cert in root_certificates)
        maximum_residual = max(cert.polynomial_residual for cert in root_certificates)
        iterations = tuple(cert.iterations for cert in root_certificates)

        return PiGConstructionResult(
            value=+final,
            level_count=level_count,
            maximum_level=level_count - 1,
            work_dps=work_dps,
            polygon_values=tuple(+v for v in polygon_values),
            extrapolated_diagonal=tuple(+v for v in diagonal),
            internal_last_diagonal_delta=+internal_delta,
            maximum_root_bracket_width=+maximum_width,
            maximum_root_polynomial_residual=+maximum_residual,
            root_iterations=iterations,
        )
