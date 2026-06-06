# Complex Torus Torsion Experiment

This note records what we learned while experimenting with the analytic map

`z ↦ (℘(z), ℘′(z))`

from the complex torus `ℂ / Λ` to the analytic Weierstrass cubic attached to
the square lattice `Λ = ℤ + ℤi`.

The corresponding short-Weierstrass companion is used only as an auxiliary
model for comparison with division polynomials. The main lesson of the
experiment is that different torsion orders require different numerical
diagnostics.

## Setup

We fixed:

- `τ = i`
- a square lattice `Λ = ℤ + ℤi`
- an analytic invariant truncation `r_inv = 16` for the main experiments
- an approximate comparison tolerance
  - absolute: `1e-4`
  - relative: `1e-2`

For some sections we also used larger truncations such as:

- `r_inv = 24`
- `r_inv = 40`
- `r_fun = 28`
- `r_fun = 40`
- `r_fun = 56`
- `r_fun = 80`

to see whether the behavior stabilized.

This turned out to matter in two different ways:

- `r_fun` controls how well the analytic evaluations `℘(z)` and `℘′(z)` are
  approximated for a fixed lattice
- `r_inv` controls how well the derived invariants `g₂`, `g₃`, and the
  resulting short-Weierstrass companion approximate the exact curve attached
  to that lattice

Equivalently:

- `r_fun` mainly improves the **point-side error**
- `r_inv` mainly improves the **curve-side error**

So a diagnostic that mixes `℘(z)` with a division polynomial built from
truncated invariants can fail for two different reasons:

- the point evaluation may still be moving
- the comparison curve may still be the wrong nearby curve

This distinction explains why different observables respond to different
truncation parameters:

- increasing `r_fun` improves the stabilization of
  `P = (℘(z), ℘′(z))`
- increasing `r_inv` improves the curve coefficients that define both
  `y² = 4x³ - g₂x - g₃` and the associated short-Weierstrass division
  polynomial `ψ_n`

So when we ask whether `ψ_n(x)` is close to zero, we are mixing:

- an approximate point `x = ℘(z)`
- with a division polynomial built from approximate invariants

and those two approximations do not improve at the same rate.

## What Worked Best For Each Torsion Order

### `n = 2`

For non-trivial `2`-torsion, the correct numerical signal is not the division
polynomial. The useful checks are:

- `|℘′(z)|`
- the residual of the cubic equation
- the distance from `x = ℘(z)` to the roots of
  `4x^3 - g₂x - g₃ = 0`

With `r_inv = 16` and `r_fun = 6, 10, 14`, we observed:

- for `(0,1;2)` and `(1,0;2)`:
  - `|℘′(z)|` decreases from about `2.37e-2` to `4.76e-3`
  - the cubic residual decreases from about `1.73` to `5.54e-2`
  - the distance to the nearest cubic root decreases from about `4.58e-3`
    to `1.47e-4`
- for `(1,1;2)`:
  - the same quantities also improve, but more slowly

Conclusion:

- `2`-torsion is visible numerically
- two of the three non-trivial classes converge much more cleanly than the
  third
- this is a better experiment than checking the even division-polynomial
  factor directly

### `n = 3`

For `3`-torsion, the best experiment is a convergence table:

- compute `℘(z)` at several `r_fun`
- compute `℘′(z)` at the same `r_fun`
- compare the cubic residual at those truncations

For representative primitive classes such as `(0,1;3)`, `(1,1;3)`, and
`(1,2;3)`, the changes

- `Δ℘`
- `Δ℘′`
- cubic residual

all decrease clearly when `r_fun` goes from `6` to `10` to `14`.

Conclusion:

- `3`-torsion already behaves convincingly as a geometric experiment
- the curve-side picture stabilizes much earlier than the division-polynomial
  criterion

### `n = 6`

For `6`-torsion, the even division-polynomial factor `ε₆(x)` is not a useful
primary diagnostic. Its values remain huge even when the torus-to-curve map is
behaving well.

The better experiment is structural:

- if `P` has order `6`, then `[2]P` should land near `3`-torsion
- if `P` has order `6`, then `[3]P` should land near `2`-torsion

We tested this by multiplying the torus representative `z` before applying the
analytic map:

- compare the image of `2z` against mapped primitive `3`-torsion
- compare the image of `3z` against mapped primitive `2`-torsion

With:

- `r_inv = 16`, `r_fun = 14`
- and also `r_inv = 24`, `r_fun = 28`

the distances were essentially zero:

- `[2]P` to `3`-torsion: at worst about `6e-14`
- `[3]P` to `2`-torsion: exactly `0` in the measured runs

Conclusion:

- the structural torsion relation is already numerically excellent
- this is a much better experiment than evaluating `ε₆(x)` directly

### `n = 7`

For `7`-torsion, the right question is whether the analytic evaluations
stabilize as `r_fun` grows.

Using `r_fun = 10, 14, 20, 28` on representative classes like:

- `(0,1;7)`
- `(0,2;7)`
- `(0,3;7)`
- `(0,4;7)`

we saw consistent decay in:

- `Δ℘`
- `Δ℘′`

For example, for `(0,1;7)`:

- `Δ℘` goes roughly `8.8e-5 → 4.9e-5 → 2.3e-5`
- `Δ℘′` goes roughly `1.2e-3 → 6.8e-4 → 3.3e-4`

Conclusion:

- the map is stabilizing
- `℘′` converges more slowly than `℘`
- the division-polynomial value `ψ₇(x)` is still numerically huge and should
  be treated only as a secondary diagnostic

There is an additional subtlety here that became clear only after pushing the
truncations much higher.

For the representative class `(0,1;7)`, with fixed `r_inv = 16`, we measured:

- `r_fun = 10`: `|ψ₇(x)| ≈ 1.2057e38`
- `r_fun = 14`: `|ψ₇(x)| ≈ 1.2171e38`
- `r_fun = 20`: `|ψ₇(x)| ≈ 1.2233e38`
- `r_fun = 28`: `|ψ₇(x)| ≈ 1.2264e38`
- `r_fun = 40`: `|ψ₇(x)| ≈ 1.2280e38`
- `r_fun = 56`: `|ψ₇(x)| ≈ 1.2288e38`
- `r_fun = 80`: `|ψ₇(x)| ≈ 1.2292e38`

So increasing `r_fun` alone does **not** drive `ψ₇(x)` toward zero. It quickly
stabilizes to a huge nonzero value.

That is not a contradiction with the exact mathematics. The exact statement is
that for the exact lattice and the exact curve,

`ψ₇(℘(z)) = 0`

for a genuine `7`-torsion class. But in the numerical experiment we are really
evaluating a division polynomial for an *approximate* short-Weierstrass
companion built from truncated invariants. If `r_fun` grows while `r_inv`
stays fixed, then `℘(z)` may become very accurate while the comparison curve
is still the wrong nearby curve. In that regime, there is no reason for the
division-polynomial value to approach zero.

We then varied `r_inv` as well, again on `(0,1;7)`, and found the stabilized
plateau moved downward:

- `r_inv = 16`: plateau near `1.229e38`
- `r_inv = 24`: plateau near `5.573e37`
- `r_inv = 40`: plateau near `2.037e37`

So for `n = 7` the two radii play clearly different roles:

- `r_fun` mainly controls stabilization of `℘(z)` and `℘′(z)`
- `r_inv` mainly controls which nearby short-Weierstrass curve is used in the
  division-polynomial comparison

In more concrete terms:

- when `r_fun` grows, the numerical image of the torsion class moves less and
  less, so the point `P = (℘(z), ℘′(z))` stabilizes
- when `r_inv` grows, the coefficients of the comparison curve move less and
  less, so the division polynomial is being built for a better approximation
  to the correct curve

This is why the experiments showed the following split behavior:

- `r_fun` made the torus-to-curve map look better
- `r_inv` made the division-polynomial value smaller

Those are not competing effects. They are corrections to two different pieces
of the same comparison.

The practical lesson is that direct `ψ₇(x)` evaluation is a badly conditioned
observable at educational truncation sizes. It is still true that the exact
limit should vanish when both approximations are taken to the exact torus and
exact curve, but it is not a realistic primary diagnostic for the finite
radii used in this experiment.

## Main Takeaways

1. The torus-to-curve map is numerically much more reliable than direct
   division-polynomial vanishing tests at moderate truncation sizes.
2. Low-order torsion should be tested using structure-adapted signals:
   - `n = 2`: `℘′(z) ≈ 0` and cubic-root proximity
   - `n = 3`: convergence of `℘`, `℘′`, and the cubic residual
   - `n = 6`: compatibility with `2`-torsion and `3`-torsion after
     multiplication
   - `n = 7`: stabilization across truncation radii
3. The division-polynomial comparison is still pedagogically useful, but for
   `n = 6` and `n = 7` it is better treated as a secondary warning signal than
   as the main correctness check.
4. For higher odd torsion such as `n = 7`, it is important to separate the
   roles of the two truncation radii:
   - `r_fun` tells us whether the analytic point evaluation is stabilizing
   - `r_inv` tells us whether the comparison curve itself is stabilizing
   Treating them as one single "accuracy knob" hides the real numerical
   behavior.
