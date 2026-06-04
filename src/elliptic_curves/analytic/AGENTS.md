# AGENTS.md for `src/elliptic_curves/analytic`

## Module mission

The `analytic` subtree is the educational complex-analytic companion to the
algebraic elliptic-curve modules.

It should grow carefully: first with small numerical domain types that state
their invariants honestly, and then with lattice sums, modular-normalization
helpers, and explanatory reports built on top of those types.

## Design rules

- Keep floating-point conventions explicit in constructors and rustdocs.
- Prefer small validated value objects for analytic knobs such as tolerances,
  truncation policies, or normalization choices.
- When a truncation object claims to represent a finite lattice-sum policy,
  make the stored parameters private and validate them at construction time.
- For square-box lattice-sum truncations, document directly that the indexing
  region is `-r ≤ m ≤ r`, `-r ≤ n ≤ r` in `ℤ²`, not a geometric disc in `ℂ`.
- When the origin matters analytically, expose the full-box and punctured-box
  counts separately instead of relying on callers to remember a hidden `-1`.
- Reject a zero truncation radius when the intended analytic routine would
  otherwise collapse to the origin-only box and stop being pedagogically
  meaningful.
- For lattice sums with mathematically meaningful convergence thresholds, use
  a dedicated typed error instead of recycling unrelated precision or
  truncation failures.
- When comparing two truncations with semantic roles like “small” and “large”,
  validate that the ordering really holds and report a dedicated typed error
  if the comparison request is malformed.
- For analytic lattice invariants, document explicitly which quantities depend
  on the scaling of `Λ` and which ones are homothety-invariant, especially
  when exposing `j` next to `g₂`, `g₃`, and `Δ`.
- When an analytic presentation and an existing algebraic model share the same
  geometric point shape, prefer a thin wrapper or type alias over duplicating
  a second point enum with the same `Infinity`/affine split.
- When approximate curve membership is exposed publicly, prefer a structured
  report helper alongside the boolean predicate so callers can inspect lhs,
  rhs, residual error, and tolerance explicitly.

## Testing expectations

- Test the positive constructor path and the typed error path for invalid
  analytic value objects.
- Keep preset constructors explicit in tests so later contributors can see
  the intended educational scale directly.
