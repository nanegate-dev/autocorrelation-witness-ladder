"""Check every certificate in this repository, from the heights up.

    python verify.py                      # all of them
    python verify.py certificates/ours_evolved_m370_n481_W1.3.json

No dependencies. Python 3.9 or newer, standard library only — `fractions` and
`json`. There is no floating point anywhere in the check: a certificate either
is a rational lower bound or it is not, and that question has an exact answer.

WHAT IS BEING CLAIMED. For a non-negative f on the line, let

    a_f(t) = integral f(x) f(x+t) dx

and let C be the smallest constant with min over t in [0,1] of a_f(t) <= C ||f||_1^2
for every such f. Because C is the SMALLEST such constant, C is the supremum of
the ratio, so EVERY f gives a LOWER bound and no f can give an upper one. A
certificate here is one f, written as a step function, together with its ratio
as an exact fraction.

WHY A STEP FUNCTION LOSES NOTHING. With f constant on each cell [i*d, (i+1)*d),
a_f is piecewise linear in t with corners only at multiples of d. The minimum
over the closed interval [0,1] is therefore attained at a corner, or at t = 1
where the value follows by exact linear interpolation between the two corners
either side. So the discretisation is not an approximation of the minimum: it
computes it.

WHAT THIS SCRIPT DOES NOT DO. It does not check that the search that produced a
witness was sensible, or that a better witness does not exist. It recomputes the
number from the heights and says whether the stored claim matches. That is the
only thing a certificate is for.

The certificates from other authors carry their provenance in the file. They are
checked by exactly the same code as ours, which is the point of keeping them
here — a verifier that only accepts its own author's files proves nothing.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: Barnard and Steinerberger, "A Lower Bound for the Autocorrelation of
#: Non-Negative Functions" — the published bound this work started from.
PUBLISHED = "0.37"

#: The upper end of the interval the constant is known to lie in:
#: 1 / (2 - 2 cos x) at the first positive root of tan x = x. Not proved here
#: and not needed for anything below; quoted so a reader can see how much of
#: the gap is left rather than only how much has been closed.
CEILING = 0.410767488189405


def ratio_of(heights: list[F], d: F) -> tuple[F, F, F]:
    """(ratio, min over [0,1] of a_f, ||f||_1) — all exact, all recomputed.

    Nothing stored in the certificate is consulted except the heights and the
    cell width. That is deliberate: a checker that reads back the answer it is
    meant to be checking is not a checker.
    """
    n = len(heights)
    # a_f at each corner k*d. The k-th is d times the correlation of the height
    # vector with itself at lag k; past the last lag it is zero.
    nodes = [
        d * sum((heights[i] * heights[i + k] for i in range(n - k)), F(0))
        for k in range(n)
    ] + [F(0)]

    inv = F(1) / d
    whole = inv.numerator // inv.denominator
    top = min(whole, n)
    candidates = [nodes[k] for k in range(top + 1)]

    # When 1/d is not a whole number, t = 1 falls between two corners and the
    # value there is the exact linear interpolation. Skipping it would report a
    # minimum that is too high, and a bound that is too high is a wrong bound.
    if inv.denominator != 1:
        if whole >= n:
            candidates.append(F(0))
        else:
            part = inv - whole
            candidates.append(nodes[whole] * (1 - part) + nodes[whole + 1] * part)

    lowest = min(candidates)
    l1 = d * sum(heights, F(0))
    if l1 == 0:
        raise ValueError("the witness is identically zero")
    return lowest / (l1 * l1), lowest, l1


def check(path: Path) -> tuple[bool, str, F]:
    cert = json.loads(path.read_text(encoding="utf-8"))
    heights = [F(x) for x in cert["heights_over_denom"]]
    d = F(cert["d"])

    if any(h < 0 for h in heights):
        return False, "a height is negative; f must be non-negative", F(0)

    ratio, lowest, l1 = ratio_of(heights, d)

    claimed = cert.get("exact_ratio")
    if claimed is not None and F(claimed) != ratio:
        return False, f"stored exact_ratio disagrees: {F(claimed)} vs recomputed", ratio
    stored_float = cert.get("float")
    if stored_float is not None and abs(float(ratio) - float(stored_float)) > 1e-12:
        return False, f"stored float {stored_float} is not this ratio", ratio
    # These two are the pieces the ratio is built from; a certificate that
    # carries them has to carry the right ones.
    if cert.get("min_g") is not None and F(cert["min_g"]) != lowest:
        return False, "stored min_g disagrees with the recomputed minimum", ratio
    if cert.get("L1") is not None and F(cert["L1"]) != l1:
        return False, "stored L1 disagrees with the recomputed norm", ratio

    return True, "", ratio


def main() -> int:
    argv = sys.argv[1:]
    paths = [Path(a) for a in argv] if argv else sorted((HERE / "certificates").glob("*.json"))
    if not paths:
        print("no certificates found")
        return 1

    print(f"the published bound this work started from: {PUBLISHED}")
    print(f"the ceiling the constant is known to sit under: {CEILING}\n")
    print(f"{'certificate':<38} {'cells':>6} {'ratio':>14}  verdict")

    rows, ok = [], True
    for path in paths:
        try:
            passed, why, ratio = check(path)
        except Exception as error:  # a malformed file is a failure, not a crash
            print(f"{path.name:<38} {'':>6} {'':>14}  UNREADABLE: {error}")
            ok = False
            continue
        n = len(json.loads(path.read_text(encoding="utf-8"))["heights_over_denom"])
        rows.append((float(ratio), path.name, n, passed, why))
        ok &= passed

    for value, name, n, passed, why in sorted(rows, reverse=True):
        verdict = "verified" if passed else f"FAILED — {why}"
        print(f"{name:<38} {n:>6} {value:>14.9f}  {verdict}")

    if rows:
        best = max(rows)
        print(f"\nthe strongest lower bound in this repository: {best[0]:.9f}")
        print(f"  from {best[1]}")
        print(f"  it leaves {CEILING - best[0]:.9f} between it and the ceiling")
    print("\nall certificates verified" if ok else "\nSOMETHING DID NOT VERIFY")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
