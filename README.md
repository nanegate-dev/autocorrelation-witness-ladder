# A certified ladder of lower bounds for the minimum-autocorrelation constant

**C ≥ 0.409821093**, as an exact fraction of two 42-digit integers, with the
witness that proves it and a checker you can run in one command.

This is [problem 6](https://google-deepmind.github.io/alphaevolve_repository_of_problems/problems/6.html)
of DeepMind's AlphaEvolve repository of problems, whose own title for it is
*"An autocorrelation problem related to difference bases"*. For a non-negative
`f` on the line, write

    a_f(t) = ∫ f(x) f(x+t) dx

and let `C` be the **smallest** constant such that

    min over t in [0,1] of a_f(t)  ≤  C ‖f‖₁²

holds for every such `f`. Because `C` is the smallest such constant, it is the
supremum of the ratio — so exhibiting any single `f` proves a **lower** bound,
and no `f` can ever prove an upper one. Everything in this repository is one
side of that: witnesses, and a way to check them.

## Check it yourself

No dependencies. Python 3.9+, standard library only.

```bash
python verify.py
```

```
certificate                             cells          ratio  verdict
ours_evolved_m370_n481_W1.3.json          481    0.409821093  verified
ours_m370_n481_W1.3.json                  481    0.409727814  verified
ours_m128_n384_W3.json                    384    0.408619504  verified
ours_m96_n288_W3.json                     288    0.407470668  verified
ours_m64_n192_W3.json                     192    0.406544583  verified
ours_m48_n144_W3.json                     144    0.405283718  verified
ours_m32_n96_W3.json                       96    0.403969677  verified
russell_n480_W1.3.json                    480    0.399213565  verified
```

`verify.py` reads nothing from a certificate but the heights and the cell
width, and recomputes the rest. There is no floating point in the check: the
ratio is built with `fractions.Fraction` and compared exactly. A checker that
reads back the answer it is meant to be checking is not a checker.

## Where the constant sits

```
0.37            Barnard & Steinerberger, published
0.399214        Russell 2026, techno-optimist/minimum-autocorrelation-bound
0.409821093     this repository
────────────────────────────────────────────────────
0.410767488     ceiling: 1/(2 − 2cos x), x the first positive root of tan x = x
```

The remaining gap is **0.000946396**.

The ceiling is quoted, not proved here. It is included so a reader can see how
much room is left rather than only how much has been closed.

## The ladder

`ladder.csv`. `f` is a step function on `pieces` cells of width `1/per_unit`.

| per_unit | cells | support | exact ratio | increment |
|---|---|---|---|---|
| 32 | 96 | 3 | 0.403969677 | — |
| 48 | 144 | 3 | 0.405283718 | +0.00131 |
| 64 | 192 | 3 | 0.406544583 | +0.00126 |
| 96 | 288 | 3 | 0.407470668 | +0.00093 |
| 128 | 384 | 3 | 0.408619504 | +0.00115 |
| 370 | 481 | 13/10 | 0.409727814 | — |
| **370** | **481** | **13/10** | **0.409821093** | — |

**The curve had not flattened where it was stopped.** The increments do not
decay. A finer rung was never run, so the best this method reaches is unknown
and **0.409821093 is not claimed to be it** — it is only the best certificate
here.

The last two rows sit at support width 13/10 rather than 3, which is the width
of the reference construction, so that pair is a direct comparison rather than
a separate experiment.

## Why a step function is exact, not an approximation

With `f` constant on each cell, `a_f` is piecewise linear in `t` with corners
only at multiples of the cell width. Its minimum over the closed interval
`[0,1]` is therefore attained at a corner, or at `t = 1`, where the value is
the exact linear interpolation between the two corners either side. The
discretisation does not approximate the minimum; it computes it.

When `1/d` is not a whole number, `t = 1` is not a corner. Skipping it reports
a minimum that is too high — and a lower bound that is too high is simply
wrong. Both checkers here handle that case.

## Checked by someone else's code

A verifier written by the same people who wrote the witness is worth something,
but not much. So:

```bash
python crosscheck.py
```

runs `third_party/verify_bound.py` — Russell's, kept byte for byte as
published, LICENSE beside it — against our witness. His script hardcodes its
input path, so `crosscheck.py` builds the layout it expects in a temporary
directory and runs it untouched. Nothing is patched or injected.

```
[verify] n=481  d=1/370  W=13/10  claimed ratio=…/… = 0.409821092522245
[verify] recomputed ratio       = …/… = 0.409821092522245
[OK]  recomputed exact ratio == pinned certificate ratio
[OK]  exact ratio exceeds published 0.37 by 0.039821
[OK]  exact ratio < proven upper bound 0.411
```

**One thing to be aware of in that output.** His script ends with a line reading
`LOWER BOUND VERIFIED: C6.6 >= 2378625/5958277 = 0.39921356...`. That fraction
is a constant baked into his script — it is *his* bound, printed regardless of
which certificate was supplied. The line that describes our witness is
`recomputed ratio`, and it reads `0.409821092522245`. We would rather point
this out than have someone find it and wonder what was being hidden.

## What is here

```
verify.py            the checker. stdlib only, exact arithmetic
crosscheck.py        runs the third-party verifier on our witness, unmodified
certificates/        eight certificates, one format, ours and the reference
ladder.csv           the table above, as data
solver/              the C++ program that produced the top rung
third_party/         verify_bound.py and its LICENSE, unmodified
```

A certificate is a JSON file holding the cell width `d`, the heights as exact
fractions, and the claimed ratio, minimum and `‖f‖₁` — also as exact fractions.
Every certificate here uses that one format, including the reference one, which
is the point of keeping it: a verifier that only accepts its own author's files
proves nothing.

## What this is not

- **Not an upper bound.** No witness can give one. The ceiling above is quoted
  from elsewhere.
- **Not a proof that 0.409821093 is the best this approach reaches.** The
  ladder was still climbing when it stopped.
- **Not a claim about `C` itself** beyond the inequality. `C` is not known.

## The witnesses are sparse, and nothing asked for that

At the top rung, 316 of 481 heights are non-zero; at the first, 28 of 96. The
zeros arrive in a few long blocks rather than scattered. That is difference-basis
structure, and it was not part of the objective, the constraints, or any
instruction — it is what the search settled on.

## Credits

The reference construction and `third_party/verify_bound.py` are Russell's,
from [techno-optimist/minimum-autocorrelation-bound](https://github.com/techno-optimist/minimum-autocorrelation-bound),
used under the licence included beside the file. The published 0.37 is Barnard
and Steinerberger, *A Lower Bound for the Autocorrelation of Non-Negative
Functions*. The problem is number 6 of the 67 in the AlphaEvolve repository of
problems, which accompanies Georgiev, Gomez-Serrano, Tao and Wagner,
*Mathematical exploration and discovery at scale* (arXiv:2511.02864).

**A note on the label `C6.6`.** It appears inside the certificate files and in
the third-party verifier, and it is a legacy name: the problems are numbered as
plain integers and there is no problem 6.6. The certificates are left byte for
byte as they were produced, because changing a certificate to tidy a label is a
habit worth not starting; nothing reads that field.

## Licence

MIT, except `third_party/`, which carries its own.
