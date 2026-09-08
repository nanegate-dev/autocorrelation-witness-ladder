#!/usr/bin/env python3
"""Independent exact-arithmetic verification of the certified lower bound

    C6.6  =  sup_{f in L^1(R), f >= 0}  min_{t in [0,1]}  a_f(t) / ||f||_1^2
          >=  2378625 / 5958277  =  0.39921356...

where a_f(t) = INT_R f(x) f(x+t) dx is the AUTOCORRELATION of f (NOT the
convolution: for the asymmetric f used here the two genuinely differ).

The witness f is a nonnegative step function on a uniform grid of width d:
    f = sum_{i=0}^{n-1} h_i * 1_{[i*d, (i+1)*d)} ,   h_i >= 0 rational.

This script reads the certificate JSON (rational heights, grid width, and the
claimed exact ratio), recomputes the ratio FROM SCRATCH in exact rational
arithmetic by a code path independent of the optimizer/certifier that produced
it, and fails loudly (nonzero exit) on any mismatch. It uses only the Python
standard library (fractions, json) -- no numpy, no floating point on the
certified path.

Why the discrete min over grid nodes is the true continuous min over [0,1]:
For a step function f on a grid of width d, the autocorrelation a_f is
piecewise LINEAR with breakpoints only at integer multiples of d (each
cell-pair overlap INT 1_{[id,(i+1)d)}(x) 1_{[jd,(j+1)d)}(x+t) dx is a triangle
in t centered at (j-i)*d of half-width d, whose slope changes only at t in
d*Z). Hence the minimum of a_f over the closed interval [0,1] is attained at a
breakpoint (a grid node k*d with k*d <= 1) or at the endpoint t=1; evaluating
those finitely many exact values yields the exact continuous minimum.
"""
from __future__ import annotations
import json
import sys
from fractions import Fraction as F
from pathlib import Path


def node_autocorr(h: list[F], d: F) -> list[F]:
    """g[k] = a_f(k*d) = d * sum_i h_i h_{i+k}, exact, for k = 0..n-1."""
    n = len(h)
    g = []
    for k in range(n):
        s = F(0)
        for i in range(n - k):
            s += h[i] * h[i + k]
        g.append(d * s)
    return g


def min_over_unit_interval(g: list[F], d: F) -> F:
    """Exact min of the piecewise-linear a_f over t in [0,1].

    Candidates: every grid node k*d with k*d <= 1 (k = 0..floor(1/d), capped at
    n where a_f = 0 past the support), plus the endpoint t = 1 obtained by exact
    linear interpolation between its surrounding nodes when 1 is not itself a
    node.
    """
    n = len(g)
    gext = g + [F(0)]                      # a_f at node k=n (t = W = n*d) is 0
    inv = F(1) / d                         # = 1/d
    kmax_node = min(inv.numerator // inv.denominator, n)   # floor(1/d), capped
    cand = [gext[k] for k in range(kmax_node + 1)]
    if inv.denominator != 1:               # t = 1 is not an integer multiple of d
        k0 = inv.numerator // inv.denominator
        if k0 >= n:
            cand.append(F(0))
        else:
            frac = inv - k0                # in (0,1)
            cand.append(gext[k0] * (1 - frac) + gext[k0 + 1] * frac)
    return min(cand)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    cert_path = root / "certs" / "certificate_n480.json"
    cert = json.loads(cert_path.read_text())

    h = [F(s) for s in cert["heights_over_denom"]]
    d = F(cert["d"])
    n = len(h)
    claim = F(cert["exact_ratio"])

    print(f"[verify] certificate: {cert_path.name}")
    print(f"[verify] n={n}  d={d}  W={n*d}  claimed ratio={claim} = {float(claim):.15f}")

    # ---- admissibility of the witness ----
    if not all(hi >= 0 for hi in h):
        print("[FAIL] f is not nonnegative (some h_i < 0)")
        return 1
    if d <= 0:
        print("[FAIL] grid width d must be positive")
        return 1
    L1 = d * sum(h)
    if L1 <= 0:
        print("[FAIL] ||f||_1 must be positive")
        return 1

    # ---- independent exact recomputation ----
    g = node_autocorr(h, d)
    min_g = min_over_unit_interval(g, d)
    ratio = min_g / (L1 * L1)

    print(f"[verify] recomputed ||f||_1 = {L1}")
    print(f"[verify] recomputed min_{{[0,1]}} a_f = {min_g} = {float(min_g):.12f}")
    print(f"[verify] recomputed ratio       = {ratio} = {float(ratio):.15f}")

    ok = True
    # 1) exact rational match against the pinned claim
    if ratio != claim:
        print(f"[FAIL] recomputed ratio {ratio} != claimed {claim}")
        ok = False
    else:
        print("[OK]  recomputed exact ratio == pinned certificate ratio")

    # 2) cross-check the pinned min_g and L1 fields if present
    if "min_g" in cert and F(cert["min_g"]) != min_g:
        print(f"[FAIL] certificate min_g {cert['min_g']} != recomputed {min_g}")
        ok = False
    if "L1" in cert and F(cert["L1"]) != L1:
        print(f"[FAIL] certificate L1 {cert['L1']} != recomputed {L1}")
        ok = False

    # 3) the theorem's displayed decimal must be a TRUE lower bound (rounded down)
    displayed = F(39921356, 10**8)         # 0.39921356, the 8-dp value printed
    if ratio < displayed:
        print(f"[FAIL] exact ratio {float(ratio)} < displayed decimal {float(displayed)} "
              f"(displayed decimal is not a valid lower bound)")
        ok = False
    else:
        print(f"[OK]  exact ratio >= displayed 0.39921356 (valid rounded-down lower bound)")

    # 4) it must strictly exceed the best published lower witness 0.37 (Barnard-Steinerberger)
    if ratio <= F(37, 100):
        print(f"[FAIL] exact ratio does not exceed the published 0.37")
        ok = False
    else:
        print(f"[OK]  exact ratio exceeds published 0.37 by {float(ratio - F(37,100)):.6f}")

    # 5) it must lie below the proven upper bound 0.411 (else a bug)
    if ratio >= F(411, 1000):
        print(f"[FAIL] exact ratio >= proven upper bound 0.411 -- impossible, bug")
        ok = False
    else:
        print(f"[OK]  exact ratio < proven upper bound 0.411")

    print()
    if ok:
        print("LOWER BOUND VERIFIED:  C6.6 >= 2378625/5958277 = 0.39921356...")
        return 0
    print("VERIFICATION FAILED.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
