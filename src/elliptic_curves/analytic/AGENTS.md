# AGENTS.md for `src/elliptic_curves/analytic`

## Module mission

The `analytic` subtree is the educational complex-analytic companion to the
algebraic elliptic-curve modules.

It should grow carefully: first with small numerical domain types that state
their invariants honestly, and then with lattice sums, modular-normalization
helpers, and explanatory reports built on top of those types.

## Design rules

- Keep floating-point conventions explicit in constructors and rustdocs.
- For Rust imports inside `analytic/`, avoid nearby relative imports by
  default. Prefer crate-root barrels such as `crate::elliptic_curves::{...}`
  when they already expose the needed surface, and otherwise use an absolute
  crate path such as `crate::elliptic_curves::analytic::...` or
  `crate::numerics::...`.
- Keep the public analytic API stratified:
  - core value objects and main algorithms may be re-exported from
    `elliptic_curves::analytic`
  - visualization stays under `crate::visualization::elliptic_curves::analytic`
  - experiment-only bundles and test-only lab reports should not leak into
    the stable public surface just for convenience
  - internal comparison payloads or helper traits should not be exposed when
    ordinary inherent methods on the public reports are enough
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
- For modular-group matrices acting on `τ`, prefer one validated value object
  for `SL_2(ℤ)` with private entries, explicit generators such as
  `S` and `T`, and checked integer arithmetic for determinant validation,
  composition, and inversion. Be explicit in docs that the stored `i128`
  model is educationally bounded, so overflow is reported honestly instead of
  being silently wrapped.
- For modular-invariance experiments such as comparing `j(τ)` against
  `j(γτ)`, prefer a structured report that stores the original point, the
  transformed point, the modular matrix, both approximations, the residual,
  the truncation, and the tolerance verdict. Document explicitly that finite
  square-box lattice truncations are coordinate-dependent, so a nonzero
  residual at fixed radius may reflect truncation error rather than genuine
  failure of modular invariance.
- For reduction to the standard fundamental domain, keep “why a modular step
  was applied” separate from “how the overall reduction attempt ended”. Use a
  step-reason enum for actual transformations and a separate terminal-status
  enum for outcomes such as already reduced or step limit reached.
- For `proptest` coverage of modular actions, prefer at least one strategy
  that samples `γ ∈ SL_2(ℤ)` directly from pseudo-random coprime integer data
  plus a Bézout witness, instead of relying only on short words in the
  generators `S` and `T`.
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
- When several analytic reports compare two complex quantities under a
  tolerance, prefer one shared composition layer such as
  `ComplexDifferenceReport` plus `ComplexApproxComparison` over repeating
  ad hoc `lhs`/`rhs`/`difference`/`tolerance` storage in each report. Keep the
  mathematically specific names (`lhs`, `rhs`, `original_j`, `transformed_j`,
  etc.) as thin accessors on top of that shared payload.
- For higher-level analytic experiment bundles that all carry the same ambient
  `τ` and lattice `Λ_τ`, prefer a tiny context trait such as
  `HasAnalyticLatticeContext` rather than forcing those reports into one
  oversized common struct.
- For high-level analytic “lab reports” that bundle several already-public
  components, add tests for internal coherence, not only for each subpiece in
  isolation. In particular, verify that duplicated views such as
  `τ ↔ Λ_τ ↔ q`, or `j` across invariant, analytic-curve, and short-model
  surfaces, stay mutually consistent inside the aggregate report.
- For Abel-Jacobi inverse-uniformization reports, prefer exposing the
  numerically meaningful decomposition of the contour integral itself:
  initial square-root branch choice, finite-segment contribution,
  compactified-ray contribution, and asymptotic tail correction. Keep
  node-by-node branch traces internal unless debugging needs force a heavier
  public report later.
- For period-recovery scaffolding, prefer reusing `ComplexLattice` for the
  recovered basis and `ComplexApproxComparison` for the recovered-`j` versus
  curve-`j` residual, instead of introducing parallel ad hoc storage for
  `ω₁`, `ω₂`, `τ`, and `close`.
- For inverse-uniformization validation reports that compare a recovered
  upper-half-plane parameter against a curve-side invariant, prefer storing
  the explicit `τ`, its standard lattice `Λ_τ`, the recovered analytic
  invariants, and the shared `ComplexApproxComparison` together so callers can
  inspect more than just the final `j` residual.
- When inverse-uniformization validation compares recovered lattice
  invariants against a target curve, keep the distinction explicit between
  direct agreement of the scale-sensitive invariants `g₂, g₃, Δ` and mere
  agreement of the modular class through `j`. If the report classifies a case
  as “same modular class but scale-sensitive mismatch”, document directly that
  this can reflect a homothety-normalization mismatch rather than a failure of
  modular recovery.
- If the Abel-Jacobi inverse layer grows beyond one medium-sized source file,
  prefer an `abel_jacobi/` module directory with focused pieces such as
  config/metadata, report types, contour selection, and integration helpers,
  while keeping the public orchestration functions in `abel_jacobi/mod.rs`.
- If a recovered period basis gets its own public wrapper type, prefer storing
  one validated `ComplexLattice` internally and deriving `ω₁`, `ω₂`, `τ`,
  oriented area, and covolume from that single source of truth. If a higher-
  level report explains how those periods were obtained from Legendre data,
  prefer wrapping the reduction report plus the integral report instead of
  duplicating raw intermediate fields without context.
- For current analytic period recovery, be explicit about the distinction between
  Legendre-side half-period integrals and the full period lattice used by
  `℘`. If `K(λ)` / `K(1-λ)` are first assembled on the Legendre side, make
  sure the public `RecoveredPeriodBasis` stores the full lattice periods, not
  semiperiods disguised as a `ComplexLattice`.
- For current analytic cubic-root recovery near the equianharmonic or otherwise
  near-pure-cubic regime, do not force the generic Cardano branch-consistency
  check when `|p|` is numerically tiny relative to the natural `|q|^{2/3}`
  scale. Prefer a documented hybrid route: pure-cubic cube-root seeds first,
  then Newton polishing.
- Keep current analytic period-recovery work under a dedicated `periods/` module
  directory. When that surface grows, prefer focused siblings such as
  recovery, normalization, reporting, or tests over re-accumulating one large
  catch-all `mod.rs`.
- Within that `periods/` subtree, keep period-lattice wrappers and
  period-side validation reports there even when they compare recovered `j`
  against curve-side invariants. A report about “did these recovered periods
  recover the right modular class?” belongs to period recovery, not to
  point-level inverse uniformization.
- When the analytic inverse direction grows beyond a single validation helper,
  prefer a dedicated sibling module such as `inverse_uniformization/` rather
  than continuing to place `τ`-validation or Abel-Jacobi point recovery under
  `periods/`. Period recovery and inverse uniformization depend on each other,
  but they are not the same mathematical stage.
- For shared period-recovery numerical knobs, prefer one validated config
  value object with private fields, explicit accessors, and educational/strict
  /loose presets over exposing a mutable bag of public counters.
- If period recovery grows a canonical modular-normalization layer, keep its
  step budget in that same validated config object instead of hardcoding an
  internal magic number inside the canonicalization helper.
- For raw complex AGM primitives, prefer a dedicated local config over reusing
  the full period-recovery config directly. If the AGM trace is exposed
  publicly, record the principal square root, the selected sign branch, and
  the resulting next-step gap explicitly so later Legendre and elliptic-
  integral layers can explain branch choices without recomputing them.
- For complete elliptic integrals in current analytic, prefer exposing both a raw
  `from_m` surface and a semantically richer `from_lambda` surface. If the
  complementary quantity is also public, keep that complement explicit in the
  function names rather than burying `1-m` or `1-λ` as a hidden convention.
- For point-level inverse uniformization in current analytic, prefer a two-level
  API: one helper that starts from an already recovered period basis and one
  end-to-end wrapper that first recovers periods from the curve. The main
  mathematical output should be a torus class in `ℂ / Λ`, so prefer returning
  a `ComplexTorusPoint` plus one chosen reduced representative rather than
  pretending the inverse lands canonically in bare `ℂ`.
- If the point-level inverse-uniformization surface grows an explicit
  Abel-Jacobi quadrature stage, prefer exposing that integral approximation as
  its own public value object before the later “reduce modulo `Λ` and validate
  against `(℘, ℘′)`” stage. That keeps the distinction visible between
  approximating the integral and interpreting it as a torus class.
- If the Abel-Jacobi implementation needs only generic complex-plane path
  geometry, such as straight segments, rays, or compactifying maps toward
  infinity, prefer placing that reusable geometry under `numerics/` rather
  than inside `analytic/inverse_uniformization/`.
- If an analytic quadrature routine uses a classical rule such as composite
  Simpson, prefer placing the generic interval-and-weight logic under
  `numerics/` and keeping only the domain-specific sampled integrand, branch
  tracking, or contour policy inside the analytic module.
- For the current Abel-Jacobi implementation, prefer transporting
  finite points first to the chosen Legendre reduction and performing the
  quadrature there, rather than integrating directly in the original
  `x`-coordinate. That keeps the branch locus visible as `{0, 1, λ, ∞}` and
  reuses the already-public affine normalization data.
- Under the current convention
  `z = ∫_x^∞ dt / sqrt(4t^3 - g₂ t - g₃)`,
  initialize the square-root branch at the finite input point using the sign
  opposite to the supplied `y`-coordinate. Document that sign choice
  explicitly whenever the code or examples discuss why the recovered value is
  `z` rather than `-z`.
- If the finite-point Abel-Jacobi quadrature uses a deterministic contour,
  prefer one simple `segment + ray` policy in the Legendre `X`-plane as the
  initial `LegendreContourStrategy` case. Keep the contour deterministic,
  explain how its angle is chosen relative to the singular locus, and expose
  both the selected strategy and the concrete contour choice in reports. At
  minimum, keep visible the start point, anchor point, chosen angle, anchor
  radius, sampled tail length, and the minimum sampled distance to the branch
  locus.
- If that contour-choice heuristic samples the segment or the compactified
  ray, do not leave those sample counts as hardcoded magic numbers. Prefer
  explicit `AbelJacobiConfig` knobs such as `segment_samples` and
  `ray_samples`, and document clearly that they tune contour scoring rather
  than the Simpson quadrature budget itself.
- When a point-level inverse-uniformization routine finishes by reducing a raw
  integral modulo the recovered lattice, prefer validating that torus
  representative by reusing `map_torus_point_to_curve(...)` and reporting the
  resulting `x`/`y` residuals. Do not treat that final validation as optional
  hidden glue.
- If callers need to compare two complex representatives only modulo an
  approximate recovered lattice, prefer a reusable helper under `lattice/`
  that searches over a small explicit box of lattice shifts and reports the
  best residual as a `ComplexApproxComparison` payload, rather than burying
  that logic inside one example.
- If that final Abel-Jacobi validation grows richer than two residual norms,
  prefer a dedicated roundtrip-validation report that keeps the recovered
  curve point, the forward-validation truncations, and the `x`/`y`
  comparison objects together instead of burying them in metadata scalars.
- If a public point-roundtrip experiment is exposed on top of the internal
  Abel-Jacobi recovery layer, prefer reusing the already recovered torus-side
  report instead of rebuilding a second parallel notion of `z_P`, contour,
  or torus class. Keep the forward-validation policy explicit as its own
  config so callers can experiment with `(wp, wp')` truncations separately
  from the inverse quadrature budget.
- Keep the forward-validation truncation policy explicit and separate from the
  inverse quadrature budget. Even for educational presets, callers should be
  able to tighten or loosen the roundtrip check without implicitly changing
  the inverse integral.
- For the current pedagogical surface, it is acceptable to reject or defer
  branch-point inputs with `y ≈ 0` instead of pretending the same quadrature
  handles them robustly. If that limitation remains, document it honestly as a
  current numerical boundary of the API.
- If that point-level inverse-uniformization surface exposes numerical
  diagnostics, keep the Abel-Jacobi metadata separate from the broader
  period-recovery metadata. Branch adjustments, lattice corrections, and
  validation residuals tell a different numerical story than AGM or Newton
  counters and should stay visible as their own typed report.
- For current analytic period-basis recovery, prefer a two-level API: one focused
  helper that starts from an already chosen Legendre reduction and one fuller
  curve-level report that also records recovered roots, the Legendre step, the
  complete-elliptic-integral report, the final basis, and the visible `τ`
  parameter.
- When both of those period-basis reports are public, prefer the fuller
  curve-level report to wrap the lower-level Legendre-to-basis report and
  expose convenience accessors, rather than duplicating the same reduction,
  integral, basis, or `τ` payload in parallel fields.
- If users frequently want only the recovered modular parameter `τ`, prefer a
  small `τ`-focused report that wraps the full period-basis recovery report,
  rather than implementing a second tau-only numerical pipeline.
- If a canonical modular representative is also exposed, prefer a third layer
  that wraps the natural `τ`-recovery report plus one explicit
  fundamental-domain reduction report. Do not silently change the meaning of
  the natural `τ` API to “already canonicalized”.
- For numerical period-recovery diagnostics, prefer one structured metadata
  value object with an explicit resolved-method enum, a status enum, and
  separate per-phase work counters over one opaque success boolean plus one
  undifferentiated iteration count. When the recovery route goes through
  Cardano branches, prefer also recording branch-selection diagnostics such as
  the Cardano discriminant, the residual of `uv ≈ -p/3`, and which branch
  indices were selected.
- For unordered complex cubic roots, do not impose an arbitrary lexicographic
  ranking just to justify `e1/e2/e3` accessors. Prefer preserving caller
  order, documenting that it has no canonical meaning, and exposing symmetric
  invariants until a later normalization introduces a mathematically
  meaningful ordering.
- For cubic-root diagnostics, keep the coarse geometric configuration
  (“three approximately real”, “one approximately real plus an approximately
  conjugate pair”, or generic complex) separate from near-collision status
  such as “nearly repeated”. Do not overload one enum to encode both ideas.
- For recovering roots of `4x^3 - g₂x - g₃`, prefer the depressed-cubic
  Cardano route with an explicit branch-consistency check `uv ≈ -p/3`, then
  Newton-polish the candidates and validate the recovered symmetric sums
  against `g₂`, `g₃`, and `e₁ + e₂ + e₃ ≈ 0`.
- When tests compare recovered cubic-root triples against expected triples,
  prefer an explicit “matches up to permutation” helper over handwritten
  chains of `||` comparisons.
- When both curve-level and invariant-level cubic-root recovery helpers exist,
  keep the invariant-level helper as the implementation surface and let the
  curve-level wrapper delegate to it.
- For cubic-root recovery reports, prefer storing `g₂`/`g₃` reconstruction
  checks as `ComplexApproxComparison` payloads rather than duplicating
  separate reconstructed-value and residual fields.
- For Legendre reduction from an unordered cubic-root triple, do not pretend
  the roots come with a canonical labeling. Prefer evaluating the six
  permutation-induced `λ` values, choosing a deterministic representative that
  stays as far as practical from the singular set `{0, 1, ∞}`, and documenting
  the tie-break rule directly in rustdocs.
- When exposing a Legendre-normalization map from
  `4(x-e₁)(x-e₂)(x-e₃)` to `X(X-1)(X-λ)`, prefer storing the affine `x`
  change-of-variables and the induced `y²` scale factor explicitly. When a
  concrete `y`-scale or invariant-differential scale is needed, choose it via
  the principal square root of `e₁ - e₂` rather than by taking a square root
  of `4(e₁-e₂)^3` directly, and document that flipping the square-root branch
  changes only a global sign.
- For Legendre reports, prefer wrapping the already computed
  `LegendreReduction` instead of duplicating roots, `λ`, scales, or
  permutations. Keep “which orbit representative was selected” separate from
  “how close the chosen `λ` is to `{0, 1, ∞}`”.
- If Legendre singularity diagnostics grow beyond one boolean, centralize the
  conditioning class, singularity-distance score, and near-zero / near-one /
  near-infinity checks in one small value object so every public helper reads
  from the same source of truth.
- For period recovery, elliptic-integral evaluation, and inverse
  uniformization, prefer dedicated typed errors such as cubic-root,
  branch-choice, validation, or Abel-Jacobi failures over collapsing distinct
  failure modes into `NumericalComparisonFailed`.

## Testing expectations

- Test the positive constructor path and the typed error path for invalid
  analytic value objects.
- Keep preset constructors explicit in tests so later contributors can see
  the intended educational scale directly.
