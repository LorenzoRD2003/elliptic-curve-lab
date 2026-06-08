# AGENTS.md for `src/visualization/elliptic_curves`

## Module mission

This subtree explains curve models, points, and group-law behavior as readable
plain text.

Its output should help a learner see both:

- the algebraic object being displayed
- the specific geometric or arithmetic rule being used

## Scope guidance

- Prefer helpers that explain short-Weierstrass curves concretely before
  generalizing.
- Keep group-law explanations explicit about cases such as identity, inverse
  points, secants, and tangents.
- For small finite groups, it is appropriate to list points and describe point
  orders, but say directly when the method is exhaustive or based on repeated
  addition.
- If a helper depends on a backend-specific capability, such as square roots
  or enumeration, surface that honestly in the explanation text.
- For short-Weierstrass isomorphism and twist explanations, keep the
  distinction explicit between:
  same `j`-invariant = isomorphic over an algebraic closure,
  versus
  isomorphic over the base field = an actual scaling witness was found in the
  current field.
- For division-polynomial helpers, keep the distinction explicit
  between:
  - rational `x`-candidates
  - lifted rational points
  - points satisfying `[n]P = O`
  - points of exact order `n`
  - comparison against exhaustive enumeration
- For Frobenius helpers, keep the distinction explicit between:
  - absolute Frobenius `π_p` metadata and relative Frobenius `π_q` metadata
  - a trace package `t = q + 1 - #E(F_q)` and the later objects derived from it,
    such as `χ_{π_q}(T)`, local zeta functions, Hasse checks, and extension counts
  - Frobenius-derived extension counts and direct exhaustive enumeration over a
    represented extension field
  - pointwise characteristic-equation terms
    `π_q(P)`, `π_q^2(P)`, `[t]π_q(P)`, `[q]P`, and the final left-hand side
  - relative-Frobenius torsion reports, where fixing is currently tautological,
    and absolute-Frobenius torsion/orbit reports, where nontrivial motion can occur
  - single-isogeny Frobenius invariance and graph-level Frobenius invariance
  - when a report contains points over extension fields, prefer their compact
    `Visualizable` / `VisualizableField` surface over raw `Debug` output, so
    the educational text still reflects quotient representatives rather than
    internal Rust structs
- For analytic helpers, keep the distinction explicit between:
  - an upper-half-plane parameter `τ` and the derived modular parameter
    `q = e^{2π i τ}`
  - one modular matrix `γ`, a transformed point `γτ`, and an accumulated
    matrix built from several reduction steps
  - finite evaluations of `℘` / `℘′`
  - the pole case at lattice points
  - the torus-side representative `z`
  - the curve-side point, membership report, and differential-equation status
  - typed torsion-vs-division-polynomial cases such as pole / odd-index /
    even-index reports
  - side-by-side modular comparisons such as `j` from Eisenstein sums versus
    `j` from a `q`-expansion
- For current analytic period-recovery helpers, keep the distinction explicit
  between:
  - one chosen period basis and a canonical lattice class
  - the naturally recovered `τ` and a later canonically reduced `τ`
  - stored cubic-root order and any mathematically meaningful classification
  - coarse geometric configuration and near-repetition diagnostics
  - a Legendre orbit label relative to input root order and any intrinsic
    data of the chosen Legendre parameter
  - reconstructed invariants `g₂, g₃` and curve-side invariants
  - inverse-uniformization validation via one explicit `τ`, the recomputed
    lattice-side invariants, and the resulting `j` residual against the
    curve-side invariant
  - whether a recovered lattice matched the curve directly at the
    scale-sensitive level `g₂, g₃, Δ` or only at the modular-invariant level
    through `j`
  - Cardano-branch selection diagnostics and final Newton/validation residuals
  - a successful-looking approximation and the numerical metadata that
    explains how it was obtained
  - a canonicalized modular representative and the accumulated modular matrix
    that produced it
  - compact pretty-printing and higher-precision diagnostic rendering when
    near-singular Legendre parameters or nearly-colliding roots would
    otherwise print misleadingly “exact” values
  - for point-level inverse uniformization, the source point `P`, the
    recovered torus representative `z_P mod Λ`, the point recovered by the
    forward map `(wp, wp')`, the forward-validation truncations, and the
    final `x` / `y` residual norms. If the public experiment reuses an
    existing Abel-Jacobi point-recovery report, keep that reuse explicit in
    the wording rather than pretending a second independent inverse path was
    recomputed.

## Formatting guidance

- Prefer compact equation strings such as `y^2 = x^3 + ax + b`.
- Use `O` for the point at infinity in compact output.
- In this subtree, prefer crate-root imports too:
  - avoid `super::...` for sibling visualization helpers unless there is a
    strong reason not to
  - prefer `crate::visualization::{...}` or other high crate-root barrels
    when they already expose the needed API surface
- In richer explanations, show the important intermediate values, such as the
  left and right sides of the curve equation or the slope used in point
  addition.
- For analytic reports, prefer explicit complex values and truncation radii
  over hiding the numerical setup behind prose alone.
- For division-polynomial explanations, include the polynomial shape explicitly:
  `polinomio en x` versus `y` times a polynomial in `x`.
- Avoid decorative ASCII art unless it materially helps understanding.

## Testing expectations

- Test the compact formatters separately from the richer descriptions.
- Test important phrases, mathematical quantities, and honest caveats.
- Include at least one finite-field example where point addition, point
  listing, and point order are all explained consistently.
- Include at least one division-polynomial example where the explanation
  mentions roots, lifted points, exact-order torsion, and the final
  enumeration comparison.
