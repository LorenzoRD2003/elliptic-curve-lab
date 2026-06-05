# Weierstrass Cubic-Root Recovery

Source: [src/elliptic_curves/analytic/periods/recovery.rs](../../src/elliptic_curves/analytic/periods/recovery.rs)

This note explains the current recovery pipeline behind
`recover_weierstrass_cubic_roots(...)` and
`recover_weierstrass_cubic_roots_from_invariants(...)`.

The goal is to recover approximate roots `e_1, e_2, e_3` such that `4x^3 - g_2 x - g_3 = 4(x - e_1)(x - e_2)(x - e_3)`.

Because the cubic is already in depressed form, the implementation uses the
classical Cardano reduction and then numerically polishes the candidates with
Newton iteration.

## High-Level Idea

Starting from the analytic curve `y^2 = 4x^3 - g_2 x - g_3`, divide the cubic equation by `4` and rewrite it as
`x^3 + px + q = 0`, where

- `p = -g_2 / 4`
- `q = -g_3 / 4`

Cardano’s ansatz writes a root as `x = u + v` with

- `u^3 = -q/2 + sqrt((q/2)^2 + (p/3)^3)`
- `v^3 = -q/2 - sqrt((q/2)^2 + (p/3)^3)`

Over `C`, each nonzero complex number has three cube roots, so the algorithm
cannot blindly take one principal branch and hope for the best. Instead, it
enumerates all three branches for `u` and all three branches for `v`, then
chooses a pair satisfying the consistency relation `uv ≈ -p/3`.

Once one consistent pair `(u, v)` is found, the three Cardano roots are

- `u + v`
- `ωu + ω^2 v`
- `ω^2 u + ωv`

with `ω = exp(2πi/3)`.

Those are then polished by Newton iteration on `f(x) = 4x^3 - g_2 x - g_3`.

## Flow Diagram

```mermaid
flowchart TB
    A["Input: analytic curve y² = 4x³ - g₂x - g₃ and recovery config"] --> B["Curve-level helper delegates to invariant-level helper with g₂ and g₃"]
    B --> C["Build validated analytic curve from the supplied invariants"]
    C --> D["Convert to depressed monic cubic x³ + px + q = 0 with p = -g₂/4 and q = -g₃/4"]
    D --> E["Compute Cardano discriminant Δ = (q/2)² + (p/3)³"]
    E --> F["Compute radicals -q/2 ± sqrt(Δ)"]
    F --> G["Enumerate the three cube-root branches of each radical"]
    G --> H["Search the 3×3 branch pairs (u, v)"]
    H --> I{"Some pair satisfies uv ≈ -p/3?"}
    I -->|"no"| Z0["Return BranchChoiceAmbiguous"]
    I -->|"yes"| J["Form the three Cardano candidates u+v, ωu+ω²v, ω²u+ωv"]
    J --> K["Newton-polish each candidate on f(x) = 4x³ - g₂x - g₃"]
    K --> L{"Residual or derivative causes failure?"}
    L -->|"yes"| Z1["Return CubicRootRecoveryFailed"]
    L -->|"no"| M["Build WeierstrassCubicRoots and validate pairwise distinctness"]
    M --> N{"Any approximately repeated pair?"}
    N -->|"yes"| Z2["Return RepeatedCubicRoot"]
    N -->|"no"| O["Check depressed-cubic relation e₁ + e₂ + e₃ ≈ 0"]
    O --> P{"Passed?"}
    P -->|"no"| Z3["Return CubicRootRecoveryFailed"]
    P -->|"yes"| Q["Reconstruct g₂ = -4(e₁e₂ + e₁e₃ + e₂e₃) and g₃ = 4e₁e₂e₃"]
    Q --> R{"Both agree with original invariants?"}
    R -->|"no"| Z4["Return CubicRootRecoveryFailed"]
    R -->|"yes"| S["Return WeierstrassCubicRoots"]
```

## Why The Validation Matters

The raw Cardano formulas produce algebraic candidates, but the implementation
still validates them numerically for three reasons.

1. Complex cube roots are branch-dependent.
   Picking inconsistent branches for `u` and `v` can break the identity
   `uv = -p/3`, so the resulting `u + v` would not actually solve the cubic.

2. Floating-point roundoff perturbs exact algebraic identities.
   Even when the correct branch pair is chosen, the initial candidates may be
   only approximate, especially after square-root and cube-root evaluation.

3. Distinct roots are part of the non-singular story.
   A repeated root would mean the cubic is colliding with a singular regime,
   so the constructor explicitly rejects approximately repeated triples.

The explicit check `e_1 + e_2 + e_3 ≈ 0` is mathematically the same as
verifying that the `x^2` coefficient vanishes, which is the signature of a
depressed cubic. It is a cheap structural sanity check that the recovered
triple still matches the original shape of `4x^3 - g_2 x - g_3`.

## Newton Polishing

Each Cardano candidate is refined with Newton iteration applied to

- `f(x) = 4x^3 - g_2 x - g_3`
- `f'(x) = 12x^2 - g_2`

The iteration stops successfully when either

- the residual is already approximately zero, or
- the Newton step becomes tiny and the post-step residual is approximately
  zero.

It fails if

- the derivative becomes approximately zero before convergence, or
- the iteration budget in `config.newton_max_iterations()` is exhausted
  without reaching an approximate root.

## Error Surface

The current recovery path can fail through these mathematically meaningful
errors:

- `BranchChoiceAmbiguous` if no branch pair satisfies `uv ≈ -p/3`
- `CubicRootRecoveryFailed` if Newton polishing or final validation fails
- `RepeatedCubicRoot` if the final triple is approximately non-distinct

The invariant-level helper may also fail earlier if `g_2, g_3` do not define
a valid non-singular analytic Weierstrass curve.

## Complexity

The implementation documents complexity as `Θ(n)` where `n = config.newton_max_iterations()`.
That estimate comes from:

- constant work to compute `p`, `q`, the discriminant, and the radicals
- constant work to search the `3 × 3` Cardano branch pairs
- at most `3n` Newton updates, one lane for each root candidate

So asymptotically the Newton polishing phase dominates the recovery routine.
