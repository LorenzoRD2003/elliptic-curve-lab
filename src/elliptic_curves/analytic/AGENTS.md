# AGENTS.md for `src/elliptic_curves/analytic`

## Module mission

The `analytic` subtree is the educational complex-analytic companion to the
algebraic elliptic-curve modules.

It should grow carefully: first with small numerical domain types that state
their invariants honestly, and then with lattice sums, modular-normalization
helpers, and explanatory reports built on top of those types.

## Design rules

- Keep floating-point conventions explicit in constructors and rustdocs.
- Once an analytic error surface becomes a stable public boundary, prefer
  implementing `Display` and `std::error::Error` with short mathematically
  honest wording.
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
- When different analytic constructions use similarly shaped truncation knobs,
  prefer distinct validated value objects when their mathematical roles may
  diverge later, even if they currently store the same primitive parameter.
- For lattice-periodic elliptic functions with poles at `Λ`, prefer reducing
  the input to the canonical fundamental-parallelogram representative first,
  then make the pole check explicit against that representative.
- When exposing first truncated evaluations of `℘`, `℘′`, or related
  elliptic functions, prefer a small structured report that stores the input
  point, the approximate value, the truncation used, and the effective term
  count.
- When several elliptic functions share the same evaluation pattern, prefer a
  generic elliptic-function trait with one default reduction-and-summation
  routine over
  copy-pasting the pole checks and lattice traversal in each function.
- Keep that shared evaluation trait internal unless users actually need to
  implement new elliptic-function families outside this module tree.
- If users do need that extension hook, prefer exposing one small public
  helper function with callbacks over promoting the whole internal evaluation
  trait to the public API.
- When multiple approximation reports store the same four core fields, prefer
  one trait helper plus default accessors over repeating identical getter
  bodies in each report impl.
- When a truncated elliptic function naturally has poles at `Λ`, prefer a
  small companion capability trait such as `HasPoleDistance` instead of
  inflating the base approximation trait for every future function.
- When mapping `ℂ / Λ` to an analytic Weierstrass curve via `(℘, ℘′)`, treat
  lattice points as the projective point at infinity instead of reporting
  them as evaluation errors.
- For reports that verify the differential equation `℘′² = 4℘³ - g₂℘ - g₃`,
  prefer reusing the existing torus-to-curve map and curve-membership report
  rather than recomputing a second inconsistent notion of lhs/rhs residual.
- When the dominant work is a full square-box traversal in the truncation
  radius `r`, prefer documenting that as `Θ(r²)` directly in rustdocs.
- For composed analytic routines that combine invariant and elliptic-function
  truncations, prefer documenting complexity as `Θ(r_inv² + r_fun²)` when
  those traversals dominate the work.
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
