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
- If a Weierstrass companion such as `ζ` is added alongside genuinely
  elliptic functions, document explicitly that it is only quasi-periodic and
  must not reuse any “reduce `z` modulo `Λ` before evaluation” shortcut for
  its actual value, even if the same reduction is still useful for pole
  detection. Prefer placing that code outside the `elliptic_functions/`
  submodule so the directory structure matches the mathematics.
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
- For modular `q`-parameters, prefer a small value object that stores both the
  validated upper-half-plane input `τ` and the derived complex number
  `q = e^{2π i τ}`. Keep pedagogical prose such as “why `|q| < 1`” in
  `visualization/`, not in the core analytic module. If that object has a
  single natural construction path from `τ`, prefer an inherent constructor
  such as `from_tau(...)` over a parallel free function.
- For early `q`-expansion experiments, prefer a dedicated validated
  `QExpansionTruncation` value object over a raw `usize`, and document
  explicitly whether `terms` counts only the nonnegative powers
  `q^0, q^1, ...` or also the principal part `q^{-1}`.
- When exposing the holomorphic Eisenstein-series family through `q`-expansions,
  prefer one validated weight object `k` with rules like “even and at least 4”
  over separate duplicated `E₄` / `E₆` family types.
- When a short explicit coefficient table is part of the educational API,
  it is acceptable to expose a small helper that returns those coefficients
  directly, as long as the docs say clearly whether the principal term is
  included or omitted.
- If a coefficient-table value object is exposed publicly, prefer storing the
  starting exponent explicitly so callers do not have to guess whether a table
  begins at `q^{-1}`, `q^0`, or another power.
- If that coefficient-table value object is also the place where exact integer
  data crosses into floating-point complex arithmetic, prefer encapsulating the
  truncated evaluation as an inherent method so that exact-to-numeric boundary
  stays localized and easy to evolve later.
- When a shared `q`-expansion abstraction covers both genuine modular forms
  such as `E₄`, `E₆` and modular functions such as `j`, prefer a neutral trait
  name like `ModularQExpansionFamily` over a mathematically narrower name like
  `ModularForm`.
- When a holomorphic Eisenstein `q`-expansion family grows beyond the special
  cases `E₄` and `E₆`, prefer one validated weight object `k` plus one runtime
  family object `E_k(q)` over duplicating separate per-weight family types.
- For that holomorphic Eisenstein family, validate directly that the weight is
  even and at least `4`, and exclude `E₂` explicitly rather than silently
  pretending it belongs to the same holomorphic modular-form surface.
- When exact `q`-series coefficients are mathematically rational in general,
  such as the Eisenstein coefficients coming from `-2k / B_k`, prefer storing
  the shared coefficient table in exact rational form instead of collapsing it
  prematurely to machine integers just because small examples like `E₄` and
  `E₆` happen to be integral.
- When a truncated modular or Eisenstein coefficient table needs every
  divisor-power sum `σ_r(n)` up to some bound `N`, prefer a shared batched
  numerics helper over recomputing each `σ_r(n)` independently inside the
  analytic module.
- Keep the exact-to-approximate boundary at coefficient-table evaluation time:
  build `q`-expansion coefficients exactly first, then convert to `Complex64`
  only when evaluating the truncated series numerically at a concrete `q`.
- If an analytic family object carries runtime parameters such as an Eisenstein
  weight, prefer verb names like `evaluate_tau(...)` over `from_tau(...)` for
  its main evaluation method so the API reads as “evaluate this family at τ”
  rather than as a constructor with hidden ambient state.
- When two analytic routes approximate the same modular quantity, prefer one
  structured comparison report that records both approximations, their
  difference, the truncations used, and the tolerance verdict, instead of
  returning only a bare boolean.
- If `q_expansion` grows beyond one tiny file, prefer a `q_expansion/` module
  directory with focused pieces such as modular-parameter types, truncations,
  specific series families, and a dedicated `tests.rs`.
- For torus-side analytic torsion, document the bridge
  `E[n] ≅ (1/n)Λ / Λ` directly in the public rustdocs so the connection to
  later algebraic `n`-torsion APIs stays visible.
- When exposing reduced torus torsion indices `(a, b; n)`, validate them at
  construction time and keep the stored fields behind accessors.
- When distinguishing primitive torus `n`-torsion, state explicitly that the
  current criterion is `gcd(a, b, n) = 1`, equivalently exact torus order `n`.
- When mapping torus torsion to the analytic cubic, document explicitly that
  the identity torsion class maps to the point at infinity, since `℘` and
  `℘′` have poles at lattice points.
- When comparing analytic torsion against division polynomials through
  `x = ℘(z)`, prefer storing an explicit even-index branch report
  (`y ≈ 0`, `ε_n(x) ≈ 0`, both, or neither) instead of only burying that
  subtlety in a warning paragraph.
- When one analytic comparison surface naturally splits into disjoint cases
  such as pole / odd-index / even-index behavior, prefer an enum with
  case-specific report structs over one catch-all struct full of `Option`
  fields.
- If the analytic torsion bridge grows across torus indices, torus-to-curve
  mapping, and division-polynomial comparison, prefer a `torsion/` module
  directory with focused subfiles plus a dedicated `tests.rs` over one large
  catch-all source file.
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
