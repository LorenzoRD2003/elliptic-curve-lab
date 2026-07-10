# AGENTS.md

@RTK.md

## Confirmed External Tools

100% confirmed in this environment on 2026-06-05:

- GP/PARI via `gp` at `/opt/homebrew/bin/gp`
- Singular at `/opt/homebrew/bin/singular`
- Octave at `/opt/homebrew/bin/octave`
- Z3 at `/opt/homebrew/bin/z3`
- CVC5 at `/opt/homebrew/bin/cvc5`
- Graphviz via `dot` at `/opt/homebrew/bin/dot`
- Gnuplot at `/opt/homebrew/bin/gnuplot`
- SageMath at `/usr/local/bin/sage`
- GAP at `/opt/homebrew/bin/gap`
- `elan-init` at `/opt/homebrew/bin/elan-init`

Important precision note:

- The GP/PARI interface is confirmed via `gp`. A standalone `pari` wrapper was
  not independently confirmed as a command in this environment.
- Separate `pari-seadata`, `pari-galdata`, `gap-system/gap`, and
  `gap-system/gap/gap` resources were not independently confirmed as standalone
  commands or paths, so do not assume those exact entry points without
  re-checking.

When a future task benefits from exact algebra, symbolic computation, SMT, or
plotting, prefer these installed tools over re-deriving everything manually.

## Project identity

`elliptic-algorithms-lab` is an educational Rust library for studying and
implementing algebraic, number-theoretic, and cryptographic building blocks.
The repository is intentionally being developed in stages:

- first, clear abstractions
- then, small correct implementations
- then, larger algorithms

This is not a race to maximize features. The main goal is to make the codebase
easy to read, easy to extend, and useful for learning.

## Primary goals

- Build a clean foundation for finite fields, polynomials, elliptic curves, and
  related algorithms.
- Favor explicit mathematical structure over magical APIs.
- Keep the code understandable for someone learning the subject and the Rust
  implementation at the same time.
- Make it easy to inspect, test, and visualize intermediate results.
- Prefer educational output surfaces such as textual explanations, operation
  tables, and polynomial formatting when they help someone understand the math.
- Prefer runnable examples when a design is easier to understand from a small
  end-to-end construction than from API signatures alone.
- When several curve-order algorithms coexist, prefer one curve-side public
  entry point with explicit method selection over scattering free-standing
  counting helpers across modules.
- If a point-count algorithm is currently used only by one concrete curve
  family, prefer keeping its executable logic with that family instead of
  extracting a generic helper too early.
- When adding a point-order recovery algorithm from a known multiple, prefer a
  step-report surface that records the prime-by-prime reductions instead of
  returning only the final order.
- If that workflow needs exact integer scaffolding such as cached prime-power
  ladders or normalized prime-power factorizations, prefer keeping that
  arithmetic source of truth under `numerics` and leaving only the genuine
  group-side peeling logic under `elliptic_curves`.
- If a finite-field order helper is really driven by the Hasse interval `H(q)`,
  prefer keeping that search/report surface in the Frobenius-side layer rather
  than presenting it as just another point-count strategy.
- If several point-order routes coexist for one curve family, prefer one
  curve-side `point_order_by(...)` entry point with an explicit order
  strategy enum, while preserving route-specific reports underneath.
- If a point-order route needs to factor an integer such as an annihilating
  multiple found in `H(q)`, prefer one canonical factoring surface under
  `numerics` instead of a second local helper in the curve layer.
- If one point-order route depends on finite-field counting/Frobenius data,
  prefer making that dependency explicit in the order-strategy payload rather
  than hiding the chosen group-order route behind an internal default.
- Support both finite and infinite base fields when the mathematics naturally
  calls for it, instead of assuming everything is cryptographic or finite from
  the start.
- Expose mathematically meaningful backend properties when they improve later
  APIs, such as whether a field family is algebraically closed.
- Keep capability boundaries explicit when only some backends honestly support
  an operation, as with square roots or future curve-side helpers.
- If exhaustive and heuristic routes coexist for one curve-side invariant such
  as point order or group exponent, prefer one public `..._by(...)` entry
  point with an explicit strategy enum and route-preserving reports instead of
  several loosely related helper functions.
- In report structs, prefer one canonical stored payload plus derived accessors
  over duplicating summary fields that can drift out of sync.
- When one route computes only a lower bound or heuristic candidate and a
  second route can later verify a related but different invariant, keep those
  stages as separate APIs and say explicitly which quantity is being certified.
- When adding a Mestre-style group-order route, keep it under the same
  curve-side `group_order_by(...)` umbrella as the other finite-field routes,
  and preserve whether the decisive uniqueness came from the original curve or
  from the quadratic twist in the returned report.
- When reviewing a new API surface, prefer lowering helper visibilities
  aggressively: conversion glue and internal wiring should be `pub(crate)` or
  narrower unless an external caller has a clear mathematical use for them.
- For educational cost/report surfaces, prefer one public value-object entry
  point such as `Type::for_kind(...)` or `Type::for_input(...)` over exposing a
  family of per-case free functions like `foo_cost()`, `bar_cost()`, and
  `baz_cost()` unless external callers genuinely need each helper as a separate
  stable symbol.
- In report structs, prefer deriving summary quantities from recorded steps
  instead of caching duplicate aggregates when one canonical history can serve
  as the source of truth.
- If a new strategy needs randomness or sampling while older strategies do
  not, prefer adding a sampler-aware sibling entry point and keeping the
  deterministic public method as a wrapper for the non-sampling routes.
- If a public strategy report already tells the right mathematical story,
  prefer upgrading its internal search engine in place, for example naive
  Hasse search to BSGS, rather than expanding the public report just to expose
  a lower-level implementation swap.
- Within `elliptic_curves::frobenius::hasse`, prefer keeping interval/bound
  objects at the top level of the `hasse` namespace and moving search engines
  under a dedicated `hasse::search` submodule; when a Hasse-side helper is
  naturally a property of a curve model, prefer exposing it as a trait or
  curve method rather than as a public free function.
- When both naive and BSGS Hasse-interval searches coexist, prefer one shared
  crate-private trait owned by `elliptic_curves::frobenius::hasse::search`
  rather than splitting that execution surface across unrelated trait files.
- Within that same `hasse::search` module, prefer keeping the crate-private
  trait surface in its own file and the concrete naive/BSGS engines in sibling
  helper files, so ownership stays local without mixing dispatch and algorithm
  internals.
- If Hasse microbenchmarks are specifically measuring search-engine choices,
  prefer colocating them under `elliptic_curves::frobenius::hasse::search`
  rather than at the broader `hasse` module level.
- Within `elliptic_curves::frobenius`, prefer grouping trace-derived value
  objects under `invariants/`, Frobenius metadata under `metadata.rs`,
  order/count/search routes under focused siblings such as `character_sum/`,
  `group_order/`, `extension_counts/`, and `hasse/`, and point-action stories
  under `orbit/` and `torsion/`, each with its own local `tests.rs` or
  sibling test files when the surface is large enough to justify dedicated
  coverage.
- When the Frobenius discriminant `Δ_π = t^2 - 4q` is later interpreted as
  quadratic-order data, prefer moving that factorization/order/report logic to
  `elliptic_curves::endomorphisms` and keeping `frobenius::FrobeniusDiscriminant`
  focused on the raw Frobenius-side datum plus elementary sign/classification
  accessors.
- Within `elliptic_curves::frobenius::orbit`, prefer keeping only genuinely
  orbit-level helpers and value objects there; if a helper takes
  `&ShortWeierstrassCurve<F>` and acts by transforming that concrete model or
  its points, prefer moving it onto `ShortWeierstrassCurve<F>` itself rather
  than exposing it as a free function.
- For high-value arithmetic algorithms such as Mestre, add property tests
  against an exhaustive baseline whenever the represented finite-field regime
  is still small enough to make that comparison honest.
- When one algorithm implementation starts mixing validation, setup, the main
  loop, and final report assembly, prefer extracting small local helpers for
  those phases before introducing a heavier state abstraction.
- If one Frobenius-side algorithm module starts owning several distinct search
  engines, prefer splitting those engines into dedicated sibling modules and
  keeping traits as thin dispatch layers.
- If one internal optimization surface has several orthogonal knobs, prefer a
  small internal config struct with independent fields over a single enum of
  mutually exclusive modes.
- For local algorithmic optimizations with measurable performance tradeoffs,
  prefer a nearby ignored microbenchmark over broad ad hoc timing elsewhere in
  the repo.
- When a finite-field group-order invariant such as parity can be detected
  through a low-degree quotient-ring computation, prefer that route over
  materializing huge dense polynomials like `x^q - x`.
- If a curve-side helper builds the defining cubic `x^3 + ax + b`, prefer
  making that conversion a method on the curve model, and keep generic modular
  polynomial exponentiation in `polynomials`.
- When a Hasse-interval BSGS search learns the parity of `#E(F_q)` in advance,
  prefer collapsing the search to one congruence class modulo `2` (equivalently,
  stepping by `[2]P`) instead of merely filtering matches after the fact.
- For small internal algorithm-config structs, prefer private fields plus
  narrow constructor/update helpers over exposing struct literals at call sites.
- For center-out BSGS traversals, prefer maintaining separate left/right
  frontiers from the center block over re-centering one giant-step state by
  repeated long jumps between alternating blocks.
- When benchmarking a heuristic justified by a distributional claim, prefer a
  deterministic corpus of varied instances that exhibits that distributional
  regime over a single fixed-instance microbenchmark.
- If a route is really a post-verification or certification step built from
  another algorithm's report, prefer keeping it as a separate API instead of
  folding it into the primary strategy enum for that invariant.
- For milestone-closing comparison examples, prefer one deterministic example
  that juxtaposes exact, heuristic, and search-based routes on the same curve
  rather than several tiny examples that each show only one route in isolation.
- If an optimization prototype benchmarks worse than the existing default,
  prefer documenting it as future work and reverting the default rather than
  leaving the slower prototype on the hot path.
- Prefer narrow capability traits for model-specific invariants when possible,
  instead of broadening foundational traits prematurely.
- When graph-like domain structures need model-specific witnesses, prefer an
  associated-type capability trait over leaking a concrete field parameter into
  otherwise generic APIs.

## Current project posture

- Educational first.
- Correctness before performance.
- Small public APIs before broad feature coverage.
- Step-by-step implementation before optimization.
- Update `AGENTS.md` on each implementation pass when a new local workflow
  rule, milestone convention, or verification constraint becomes relevant.
- At the start of each new conversation in this repo, before implementation or
  deep code analysis, first check for an existing Repomix artifact such as
  `repomix-output.xml` or `repomix-output.md`. If present, read it and use its
  token summary to identify the most expensive files, especially `AGENTS.md`,
  long docs, and very large source files; if absent, generate a compressed
  Repomix artifact and perform that same analysis before proceeding.
- After that Repomix check, explicitly prune context to the active work area.
  For localized work, avoid loading unrelated heavyweight guidance files,
  oversized docs, and large source files unless they are materially relevant
  to the current task. In particular, keep only the root `AGENTS.md`, the
  nearest module-local `AGENTS.md`, the active feature plan, and adjacent code
  in immediate focus unless the task actually depends on broader context.
- For Rust test verification during focused local work, prefer
  `cargo test -q` with the narrowest module-relevant filter that honestly
  covers the touched code before escalating to broader suites. Do not run the
  full `cargo test -q` for a localized edit unless the touched surface really
  crosses enough module boundaries that narrower verification would be
  misleading.
- When introducing native coordinate formulas, especially affine/projective
  group-law formulas, document the exact mathematical identities in rustdocs
  near the implementing helper or method, state the coordinate chart or model
  convention being used, and mention asymptotic complexity when that helps the
  educational story.
- Once a curve family has a native projective group-law engine, prefer making
  the affine public group-law wrappers delegate to that single executable core
  instead of maintaining a second affine execution path in parallel.
- When a staged curve family gains a reduction or transport layer, keep the
  reduction object explicit and preserve the coordinate-change data as
  first-class state instead of returning only the reduced companion model.
- For Montgomery curves `B y^2 = x^3 + A x^2 + x`, treat the model itself as
  valid in every characteristic different from `2`; any extra restriction
  involving characteristic `3` belongs to a specific conversion route such as
  Montgomery-to-short-Weierstrass normalization, not to the Montgomery model
  constructor.
- For Montgomery invariants in the repo's `c4,c6,Δ,j` normalization, keep the
  formulas consistent with `j = c4^3 / Δ` and
  `c4^3 - c6^2 = 1728Δ`; in particular use
  `c4 = 16(A^2 - 3)/B^2`, `c6 = 32A(9 - 2A^2)/B^3`, and
  `Δ = 16(A^2 - 4)/B^6`, not a lower power of `B` in the discriminant.
- For Montgomery `LiftXCoordinate` work in odd characteristic, treat the fiber
  above `x` as the square-root problem
  `y^2 = (x^3 + A x^2 + x)/B`; if `B != 0` has already been validated by the
  constructor, prefer reusing that scaled square-root story directly instead of
  routing early Montgomery lifting through another curve model.
- For the classical Montgomery-to-short-Weierstrass reduction in
  characteristic different from `2` and `3`, use the affine change of
  variables `x = B X - A/3`, `y = B Y`, equivalently
  `X = (x + A/3)/B`, `Y = y/B`, and keep the characteristic-`3` restriction
  attached to this conversion route rather than to the Montgomery model
  itself.
- For whole-curve Montgomery-to-general conversion, prefer the direct affine
  rescaling `x = B X`, `y = B Y`, which yields the general-Weierstrass model
  `Y^2 = X^3 + (A/B)X^2 + (1/B^2)X`; this route stays available in
  characteristic `3` because it does not divide by `3`.
- For whole-curve short/general-to-Montgomery conversion over the same base
  field, treat compatibility as a certified extra condition: the current
  staged route may search for a rational `2`-torsion root `α` of the short
  cubic together with a square root of `3α^2 + a`, and should fail honestly
  when no such Montgomery presentation is certified over the represented base
  field.
- For staged native Montgomery affine group-law work on
  `B y^2 = x^3 + A x^2 + x`, keep the formulas explicit in docs and tests:
  `-(x,y) = (x,-y)`,
  `λ_add = (y2-y1)/(x2-x1)`,
  `x3 = B λ^2 - A - x1 - x2`,
  `y3 = λ(x1-x3) - y1`,
  `λ_double = (3x^2 + 2Ax + 1)/(2By)`,
  `x([2]P) = B λ^2 - A - 2x`,
  `y([2]P) = λ(x-x([2]P)) - y`.
- For staged finite-field Montgomery wrappers, prefer the same honesty pattern
  as the general-Weierstrass family: exhaustive counting, point orders, Hasse
  search, and exponent accumulation may be native once `LiftXCoordinate` and
  the affine group law exist, while deeper finite-field routes such as
  quadratic-character counting or Schoof should delegate to the short
  companion only when the Montgomery-to-short reduction is actually available.
- For staged `TwistedEdwardsCurve<F>` work in characteristic different from
  `2`, prefer making the Montgomery family the canonical coefficient-level
  bridge: whole-curve conversion should be owned directly by the
  Twisted-Edwards/Montgomery pair, while short/general Weierstrass reuse should
  come from composition rather than from a second independent reduction story.
- For the first Twisted-Edwards descriptor milestone, keep the public surface
  narrow: validated coefficients, equation formatting, and classical
  invariants are enough. If those invariants are implemented through the
  Montgomery normalization formulas, leave that derivation legible in code and
  rustdocs.
- When that same Twisted-Edwards bridge is introduced, treat point transport
  honestly as a birational chart issue rather than assuming the current
  total-point `CurveModelConversion` contract applies automatically on affine
  points; resolve exceptional-point semantics explicitly before promising
  roundtrip point transport across every affine point of both models.
- If a new curve family has a finite affine identity, such as the Edwards
  neutral element `(0, 1)`, do not let shared enumeration helpers silently
  duplicate it. Prefer making the shared `EnumerableCurveModel` path
  identity-aware rather than forcing that model back into an artificial
  point-at-infinity convention.
- For the first Twisted-Edwards point-transport API, prefer names that make
  the birational-open-subset semantics explicit, such as
  `try_point_to_montgomery_open(...)`, rather than names that read like a total
  global equivalence of affine point sets.
- Keep the three Twisted-Edwards transport layers explicitly separated in docs
  and API reviews: whole-curve conversion, birational point transport on an
  affine open, and any later total rational-point correspondence.
- For the first native affine Twisted-Edwards group law on
  `a x^2 + y^2 = 1 + d x^2 y^2`, prefer the generic coefficient family with
  denominator-failure handled honestly; do not call the formulas "complete"
  unless a later restricted subfamily and proof obligations are documented
  explicitly.
- For the first Twisted-Edwards projective milestone, prefer a dedicated
  `twisted_edwards/projective/` module with an
  `ExtendedTwistedEdwardsPoint<F>` type in coordinates `(X:Y:Z:T)`
  representing `x = X/Z`, `y = Y/Z`, and `T = XY/Z`. Validate membership
  with the pair of equations `aX^2 + Y^2 = Z^2 + dT^2` and `XY = ZT`, keep
  the neutral element as `(0:1:1:0)`, and treat `Z = 0` honestly as a
  non-affine projective case rather than silently coercing it back into the
  affine identity convention.
- For Montgomery educational examples and visualization, prefer showing three
  surfaces side by side when they are available: the native Montgomery model,
  the explicit short companion, and the direct general-Weierstrass embedding,
  together with one concrete point transport or group-law comparison so the
  user can see both what changes and what stays invariant.
- Even when one direction of a model conversion is mathematically trivial,
  prefer exposing that inclusion explicitly in the same reduction layer so the
  API stays symmetric and testable.
- When multiple curve families are expected to coexist, prefer a public shared
  conversion trait plus a private witness type over exposing one family-specific
  reduction struct as part of the stable API.
- If a new shared trait naturally depends on one minimal capability such as
  `CurveModel`, it is acceptable to implement that minimal capability early for
  a staged model family without pulling in the rest of its later milestone work.
- When two error types are genuinely coupled across one abstraction boundary,
  prefer total `From` conversions over repeated `map_err` matches, as long as
  any lossy degradation is explicit and mathematically honest.
- When a whole curve can be reinterpreted as another model without extra
  caller-supplied data, prefer exposing that ergonomic curve-level conversion
  through `From` or `TryFrom` in addition to any richer witness object used
  for point transport.
- The same preference applies to tightly coupled local helper errors: if one
  model-specific helper error always degrades to one shared domain error, use a
  direct `From<LocalError> for SharedError` instead of keeping a free-standing
  mapper function at the call site.
- Treat `LiftXCoordinate` as the affine-fiber story for the projection
  `x : E -> A^1`, not as a synonym for “compute a square root of `rhs(x)`”.
  If one model recovers that fiber through square roots and another through a
  shifted quadratic or characteristic-`2` Artin-Schreier solve, keep the trait
  generic enough to describe the fiber and keep the solving route model-specific.
- For general-Weierstrass x-lifting prep, prefer one dedicated helper for the
  `y^2 + uy = v` equation over widening `SqrtField` or hardcoding the
  coordinate algebra directly into the eventual `LiftXCoordinate` impl.
- When that general-Weierstrass helper eventually feeds `LiftXCoordinate`,
  prefer one unified fiber solver that dispatches honestly between: completing
  the square in odd characteristic, inverse-Frobenius square roots when
  `u = 0` in characteristic `2`, and Artin-Schreier solving when `u != 0`.
- In that same general-Weierstrass lifting layer, prefer fiber-oriented names
  such as `y_fiber_equation`, `linear_coefficient`, and `right_hand_side`
  over operational names like “quadratic for x” when the code is really
  describing the fiber of `x : E -> A^1`.
- If one helper module in that lifting layer starts mixing equation data,
  backend solver dispatch, and curve-extension methods, prefer promoting it to
  a small submodule directory with one file per responsibility instead of
  keeping one long mixed helper file.
- Once such a helper becomes a submodule directory, prefer colocating its
  narrowly focused tests under that directory as well, instead of letting the
  parent module's catch-all `tests.rs` keep growing.
- More generally, if a curve-trait file accumulates exact integer helpers that
  do not depend on curve semantics, prefer moving them into `numerics/` and
  keeping only genuinely curve-domain validation or error helpers local.
- For staged general-Weierstrass group-law support, prefer one native affine
  negation formula `-(x, y) = (x, -y - a1*x - a3)` together with honest affine
  addition/doubling formulas before introducing projective machinery.
- Once such a family gains a staged projective capability, prefer keeping the
  affine route as an explicit oracle/bridge until the projective formulas are
  independently validated, rather than switching the public group law and the
  new projective layer simultaneously.
- Treat that affine group law as transitional. Leave an explicit TODO near the
  implementation that the long-term replacement should be a projective-coordinate
  general-Weierstrass law, even after affine formulas land.
- When a new curve family gains native group-law support, add both exhaustive
  tiny-field checks in sensitive characteristics and `proptest` coverage over
  representative supported characteristics; if the generators look reusable,
  prefer placing them under `proptest_support::elliptic_curves`.
- Once such a family also satisfies the blanket bounds for
  `EnumerableCurveModel`, `FiniteGroupCurveModel`, `FrobeniusTraceCurveModel`,
  or `IsogenyKernel::cyclic`, certify that compatibility explicitly with
  targeted tests rather than treating the blanket impls as “obviously fine”.
- For nontrivial curve families, prefer a `tests/` directory with files split
  by testing intent such as construction, reduction, point lifting, group law,
  and compatibility, instead of one growing catch-all `tests.rs`.
- When adding high-value wrapper methods on a staged curve family, document at
  each wrapper whether the route is native to that model or currently
  delegated through a companion model. Prefer native exhaustive/small-group
  routes when they already exist, and delegate only the genuinely model-specific
  algorithms that have not been generalized yet.
- When a staged curve family coexists with a companion model, close each
  compatibility milestone with a dedicated `tests/compatibility.rs` suite that
  checks invariant preservation, point transport, and operation/order
  compatibility explicitly instead of leaving those guarantees scattered across
  unrelated test files.
- When that same staged family reaches the educational/examples milestone,
  prefer a model-specific visualization helper plus one runnable example that
  shows the source equation, the companion reduction when available, and one
  concrete calculation compared across both models.
- Once that unified solver exists, let `Q` and `ComplexApprox` reuse the
  odd-characteristic path immediately, and document separately why current
  `Q`-extension backends still sit outside the trait until they gain honest
  square-root support.
- Textual explanations and visualizations are welcome when they improve
  understanding.
- The small finite-field graph layer starts with a deliberately small scaffold for educational
  `ℓ`-isogeny graphs over tiny prime fields before any broader graph or
  volcano machinery is attempted.
- For that graph layer, prefer small dedicated helper modules when one narrow piece
  of graph logic grows its own vocabulary or tests, rather than keeping every
  helper inside one builder file.
- The division-polynomial layer adds educational division-polynomial and rational-torsion
  tooling for short-Weierstrass curves over small enumerable fields.
- For that division-polynomial layer, prefer separating:
  - generic torsion-order logic under `elliptic_curves`
  - division-polynomial shape and evaluation logic under
    `elliptic_curves::short_weierstrass::division_polynomials`
  - isogeny-graph kernel wrappers under `isogenies::graphs`
  - visualization and example walkthroughs under
    `visualization::elliptic_curves` and `examples/`
- Analytic examples should present the analytic story in the same
  educational style: state `τ`, the lattice basis, the truncations chosen,
  the computed invariants, and any numerical caveats or approximation
  tolerances explicitly.
- Legendre-reduction examples should show the source cubic roots,
  the chosen Legendre parameter `λ`, the relevant `S₃` orbit/permutation
  information, and at least one explicit numerical sanity check for the affine
  normalization and scale identities.
- When a current analytic example is specifically about Legendre reduction,
  prefer also showing one “what changes / what does not change” comparison
  under root permutation and one controlled rejection for repeated-root or
  singular input, so the numerical API boundaries are pedagogically visible.
- Period-recovery examples should present the story in the same
  educational style: state the recovery config, the source or reconstructed
  invariants, the recovered roots, any root-classification report, and the
  caveat that stored complex-root order is not canonical. When printing the
  analytic cubic itself, prefer the specialized formatter under
  `visualization::elliptic_curves` over a raw generic `Display` if the latter
  keeps avoidable `+ 0i` noise or misleading sign surfaces. Prefer including
  at least one numerically harder case, for example a noisy-invariants case or
  one where Newton polishing is genuinely exercised.
- For current analytic cubic-root recovery, treat the near-pure-cubic regime
  honestly: when the depressed coefficient `p` is numerically tiny relative to
  the natural `|q|^{2/3}` scale, prefer a documented hybrid strategy that
  seeds Newton from the cube roots of `-q` instead of insisting on a fragile
  generic Cardano branch match.
- For current analytic period-basis recovery specifically, keep the distinction
  explicit between Legendre-side half-period integrals and the final full
  period lattice used by `℘`. Public period bases and end-to-end examples
  should expose full lattice periods, not semiperiods.
- When a current analytic example is explicitly end-to-end from a curve, prefer
  showing both the full period-basis recovery report and any τ-focused wrapper
  report side by side, so it stays clear that the user-facing “just give me τ”
  surface reuses the same underlying recovery pipeline.
- If that same example also shows a canonical modular representative, prefer
  showing the natural recovered `τ`, the canonically reduced `τ`, and the
  modular matrix relating them, so the distinction between recovery and
  normalization remains visible.
- Current analytic inverse-uniformization examples should present the point-level
  story in the same educational style: state the source curve point, the
  recovered period basis used as ambient lattice data, the Abel-Jacobi config,
  the chosen contour convention, including its sampled tail length when that
  convention is reported numerically, the recovered torus representative, and
  the final validation residuals after mapping back through `(℘, ℘′)`.
- When a current analytic point-roundtrip example is built from synthetic torus
  samples `z -> P`, prefer also printing whether the recovered torus class
  matches the original source class in `C / Λ`. Include at least one
  deliberately validation-limited case so the example distinguishes “inverse
  recovery looks reasonable” from “the chosen forward `(℘, ℘′)` validation
  budget was too weak”.
- In those point-roundtrip examples, prefer using a pedagogical validation
  tolerance for the “healthy” cases so the output reflects the intended
  numerical story, while reserving at least one case where the same tolerance
  fails because the forward validation budget or the geometry is genuinely
  harder.
- For the current Abel-Jacobi implementation, the intended
  pedagogical convention is:
  - transport first to Legendre form
  - integrate along one deterministic `segment + ray` contour in the
    normalized `X`-plane
  - follow the square-root branch by continuity
  - use the sign opposite to the input `y`-coordinate at the starting point so
    the convention `z = ∫_x^∞ dt / sqrt(4t^3 - g₂ t - g₃)` matches the local
    uniformization parameter
    Keep that convention explicit in docs and examples rather than burying it in
    implementation details.
- If that contour-selection layer uses diagnostic sampling of the finite
  segment or compactified ray, prefer explicit `AbelJacobiConfig` knobs over
  hardcoded constants, and explain in docs/examples that those knobs tune the
  contour-selection heuristic rather than the Simpson integration budget.
- When current analytic inverse uniformization validates a recovered torus
  representative, prefer reusing the existing torus-to-curve map and reporting
  the `x`/`y` residuals explicitly. If a finite-point example only stabilizes
  under a stricter or more explicit `AbelJacobiConfig` than the loose preset,
  show that choice honestly instead of pretending every preset is equally
  robust on every input.
- For current analytic AGM work, keep the raw complex AGM primitive separate from
  the higher-level complete-elliptic-integral API. Prefer one dedicated AGM
  config that can be derived from `PeriodRecoveryConfig`, and when exposing an
  educational trace, record the square-root branch chosen at each iteration.

The most mature parts of the repository now include `fields`, `polynomials`,
and the finite-field / analytic elliptic-curve layers, especially:

- `Fp<M, LIMBS>` and `FpElem<M, LIMBS>` for exact static Montgomery
  prime-field arithmetic
- `Q` for exact rational arithmetic over `BigRational`
- `ComplexApprox` for approximate numerical experiments over `C`
- `SqrtField` as a small capability trait for backends that can produce square
  roots honestly, with current implementations for `Fp<M, LIMBS>`,
  `Fp<M, LIMBS>`, `Q`, and
  `ComplexApprox`
  - a quadratic-character capability for finite fields of odd characteristic,
    with explicit `0 / residue / non-residue` values rather than a hidden
    ad hoc integer convention
  - a quadratic-character group-order route
    `#E(F_q) = q + 1 + \sum_x χ(x^3 + Ax + B)` kept explicit as distinct from
    exhaustive affine-point enumeration
- `PthRootExtraction` as a characteristic-`p` capability trait on values:
  finite-field elements use inverse Frobenius to expose canonical `p`-th
  roots, while denser objects such as polynomials can implement the same trait
  with stricter existence criteria
- `EnumerableFiniteField` for small finite backends that can honestly expose
  their full element set
- `ExtensionField<S>` / `ExtensionFieldSpec` as a type-level quotient-field
  design for algebraic extensions and towers over arbitrary base fields,
  including working quotient arithmetic and inversion
- `PolynomialFieldElement<F>` as an autocontained quotient-value layer with
  canonical reduction, quotient-class equality, and basic arithmetic over a
  stored modulus
- `RationalFunction<F>` as a first autocontained univariate rational-function
  value layer over `DensePolynomial<F>`, with gcd-based reduction,
  denominator-monic normalization, arithmetic, inversion, and formal
  differentiation
- `RationalFunctionField<F>` as the zero-sized field family whose element type
  is `RationalFunction<F>`, separating field metadata from stored values in the
  same style as other field-family backends
- `PolynomialModulus<F>::check_field_modulus_requirements()` as the bridge
  from polynomial irreducibility results into field-domain quotient checks
- dense, sparse, and multivariate polynomial representations over fields
- dense Euclidean division and dense gcd over fields
- formal univariate derivatives for dense and sparse polynomials through the
  shared `UnivariatePolynomial` surface
- shared univariate `gcd` support across dense and sparse polynomials, with
  the current sparse implementation delegating to the dense Euclidean backend
- baseline irreducibility classification over prime fields, plus
  field-theoretic reducibility classification for algebraically closed
  backends such as `ComplexApprox`, plus an exact but partial backend for `Q`
  that certifies some cases and returns an honest inconclusive error otherwise
- univariate evaluation plus baseline interpolation through the classical
  Lagrange formula
- a typed `PolynomialError` surface shared by polynomial-domain APIs and
  explanation helpers
- text-based visualization helpers for prime fields, rationals, polynomials,
  complex numbers, square-root behavior, and short-Weierstrass curve helpers
  ranging from point membership and addition explanations up through compact
  and verbose finite-group summaries for small enumerated curves
- a substantial `elliptic_curves` layer, centered on affine points,
  short-Weierstrass curves, discriminants, curve-membership checks,
  `x`-coordinate lifting, small-field point enumeration, explicit additive
  group-law traits and model capabilities, small-group helpers such as torsion
  checks, point orders, group order, and group exponent, classical
  short-Weierstrass invariants such as `c4`, `c6`, and `j`, plus an explicit
  Frobenius layer that distinguishes:
  - absolute Frobenius `π_p` on coordinates and Frobenius twists
  - relative Frobenius `π_q` as an endomorphism over the chosen finite base field
  - Frobenius trace data recovered from the counting formula
    `t = q + 1 - #E(F_q)` over small finite base fields
  - exact Hasse-bound verification derived from that trace through the
    equivalent integer inequality `t^2 <= 4q`
  - exact discrete Hasse intervals `H(q)` derived from the same finite-field
    order, with integer endpoint arithmetic suitable for order searches and
    Mestre-style uniqueness checks
  - ordinary versus supersingular classification derived from the same
    trace through the divisibility criterion `p | t`
  - a first relative-Frobenius-on-torsion report over exact rational
    `n`-torsion, included as a pedagogical bridge to later nontrivial
    extension-field Frobenius/torsion behavior
  - a first nontrivial absolute-Frobenius-on-torsion report for
    base-defined short-Weierstrass curves viewed over extension fields
  - explicit Frobenius orbit helpers, including nontrivial absolute orbits
    over extension fields and the trivial relative-orbit counterpart
  - fixed-point behavior that begins to separate `E(F_p)` from larger
    finite-field point sets such as `E(F_{p^r})`
  - associated metadata constructors such as
    `AbsoluteFrobenius::for_field::<F>(...)` and
    `RelativeFrobenius::for_field::<F>(...)` rather than free-standing
    module-level builders
  - an implemented educational `elliptic_curves::endomorphisms` layer derived
    from finite-field Frobenius data such as `q`, `t`, `χ_{π_q}(T)`, and
    `t^2 - 4q`, with candidate-order and quadratic-order scaffolding kept
    separate from stronger geometric claims
  - a first function-field layer for short-Weierstrass curves, modeling
    `F(E) = F(x) ⊕ yF(x)` through pairs of rational functions
    `(A(x), B(x))` representing `A(x) + yB(x)`, with multiplication reduced by
    the specific short-Weierstrass relation `y^2 = x^3 + ax + b`, plus public
    substitution helpers for evaluating polynomials and rational functions in
    the distinguished `x`-coordinate at a function-field element
  - within `elliptic_curves::short_weierstrass::function_fields`, prefer keeping methods whose
    receiver is `ShortWeierstrassFunctionField<F>` under
    `function_fields::field` only when they are genuinely ambient-field
    responsibilities such as embeddings and generic-point helpers; keep
    function-field point validation, conversion, group-law adapters, and
    point-arithmetic wrappers together under `function_fields::point`, even
    when some wrappers still take `&ShortWeierstrassFunctionField<F>` as the
    receiver
  - within `elliptic_curves::short_weierstrass::function_fields::value`, prefer keeping
    `arithmetic.rs` focused on the public algebraic operations of
    `ShortWeierstrassFunction<F>` and moving curve-compatibility checks,
    right-hand-side builders, and error-mapping helpers into a separate
    crate-private internal helper file
  - when a helper in short-Weierstrass function fields is really describing the short-Weierstrass
    cubic `x^3 + ax + b` rather than the value layer itself, prefer placing it
    on `ShortWeierstrassCurve<F>` as a crate-private helper instead of keeping
    it under `function_fields::value`
- the first usable pieces of division-polynomial torsion tooling, including:
  - generic exact-order helpers such as `point_has_exact_order(...)` and
    `points_of_exact_order(...)`
  - educational division-polynomial shape tracking through
    `DivisionPolynomialForm<F>`
  - low-degree base division polynomials `ψ_0` through `ψ_4`
  - recursive odd/even division-polynomial construction over small fields
  - pointwise and `x`-coordinate evaluation helpers
  - rational `x`-candidate, torsion-candidate, torsion-point, and exact-order
    torsion-point recovery surfaces derived from division polynomials
  - comparison reports between division-polynomial recovery and exhaustive
    torsion enumeration
- the first usable pieces of `elliptic_curves::isomorphisms`, including a
  small `CurveIsomorphism` trait plus explicit short-Weierstrass base-field
  scaling isomorphisms with cached codomains and exhaustive witness search
- substantial pieces of `isogenies`, including explicit finite kernels, Vélu
  isogenies on short-Weierstrass curves, exhaustive structural verification
  helpers, strict and bridged composition on small finite curves,
  scalar-multiplication isogenies `[n]`, exhaustive map-equality helpers, dual
  Vélu search by enumerating tiny rational kernels and testing both duality
  relations on rational points, public helpers for checking
  `\hat{\phi} \circ \phi = [deg(\phi)]` and
  `\phi \circ \hat{\phi} = [deg(\phi)]` exhaustively, graph
  construction/verification over tiny prime fields, Frobenius-relation reports
  for isogenies and graphs, together with a short-Weierstrass pullback layer on
  function fields that represents
  `\phi^* : F(E') -> F(E)` through the images of `x'` and `y'`, supports
  substitution of rational functions and `A(x') + y'B(x')`, and composes those
  pullbacks contravariantly; current short-Weierstrass Vélu isogenies can also
  export these pullbacks explicitly via `x_pullback`, `y_pullback`, and
  `as_function_field_map`; the first separability-side surface is now a
  `DifferentialPullbackReport` driven by the invariant differential multiplier;
  explicit absolute and relative Frobenius isogenies are now also modeled as
  purely inseparable maps with point evaluation, pullback formulas on function
  fields, and separate separable / inseparable degree metadata;
  `proptest_support` now also includes valid generators for these pullback maps
  and for composable pullback-map pairs
- a substantial complex-analytic layer under `elliptic_curves::analytic`,
  including validated upper-half-plane and lattice data, truncated Eisenstein /
  `℘` / `ζ` evaluation, modular actions and `q`-parameters, period recovery,
  inverse uniformization, torus torsion comparisons, and related
  visualization/report surfaces
- text-based visualization helpers for dual-isogeny workflows,
  including composition summaries, scalar-multiplication summaries, dual
  isogeny summaries, and exhaustive dual-verification reports suitable for the
  final dual-isogeny example
- runnable educational examples under `examples/`, including extension towers
  plus walkthroughs for Frobenius and group-order algorithms, group structure,
  isomorphisms, Vélu isogenies, dual isogenies, `ℓ`-isogeny-graph exploration,
  division-polynomial torsion recovery, and the analytic period-recovery /
  inverse-uniformization / torsion-comparison story that show how the APIs and
  visualization surfaces are meant to be used

## Code style expectations

- Prefer idiomatic, readable Rust over clever or excessively generic code.
- Keep modules small and focused.
- Prefer explicit naming over short cryptic names, especially in educational
  code.
- For Rust imports, prefer crate-root barrels and absolute paths over
  relative-module imports:
  - avoid `super::...` / `super::super::...` imports by default
  - prefer `crate::elliptic_curves::{...}`, `crate::visualization::{...}`,
    `crate::fields::{...}`, or another high crate-root barrel when one exists
  - when no suitable barrel exists, still prefer an absolute crate path such
    as `crate::elliptic_curves::analytic::...` over a relative import
- Public APIs should be conservative and easy to explain.
- Use `Result` for recoverable validation and arithmetic errors.
- Prefer typed domain errors such as `FieldError`, `PolynomialError`, and
  `CurveError` over raw string errors once a module has more than one
  meaningful failure mode.
- Add `///` rustdocs to public traits, structs, functions, and any non-obvious
  internal helper that carries important meaning.
- Use `todo!()` only when deferral is intentional and the message explains what
  remains undecided or unimplemented.

## Educational writing rules

- Document mathematical assumptions directly in rustdocs or nearby comments.
- For algorithmic public APIs, document asymptotic complexity when it adds
  real educational value. Prefer `Θ(...)` notation over `O(...)` for the
  current convention, and say clearly whether the estimate is counting group
  operations, field operations, full enumeration passes, memory, or
  bit-complexity.
- Explain why a representation was chosen when the choice is not obvious.
- Prefer examples and concrete terminology such as `GF(17)` or `F[x]/(m(x))`
  over abstract wording when possible.
- Avoid hiding domain invariants in “smart” helper layers; make them visible in
  types, constructors, or docs.
- When writing LaTeX in Markdown documentation or explanations, use single
  dollars `$...$` for inline math and double dollars `$$...$$` for display
  math. Do not use backticks such as `...` or fenced blocks like
  `math` for mathematical notation.
- If an implementation is approximate, pedagogical, incomplete, or not suitable
  for production cryptography, say so explicitly.

## Architecture conventions

- Keep domain boundaries clear:
  - `fields`: field abstractions and implementations
  - `polynomials`: polynomial representations and later polynomial algorithms
- `visualization`: educational text-formatting and explanation helpers split
  by mathematical domain, including both compact and verbose elliptic-curve
  group-reporting surfaces when the group is small enough to enumerate
  - `elliptic_curves`: curve models and point representations
  - `elliptic_curves::analytic::periods`: recovery of roots, Legendre data,
    complete elliptic integrals, and period bases
  - `elliptic_curves::analytic::inverse_uniformization`: validation of
    recovered `τ`/lattice data plus point-level Abel-Jacobi recovery back to
    torus classes
  - `algorithms`: reusable algorithmic building blocks
  - `utils`: project-wide helpers that do not belong to a narrower domain
- Re-export only stable, intentional entry points from `lib.rs` and `mod.rs`.
- Keep crate-root and high-level domain barrels focused on core mathematical
  and algorithmic surfaces. Do not re-export visualization helpers, ad hoc
  experiment reports, or internal comparison payloads from `lib.rs` just for
  convenience; those should stay under their natural namespaces such as
  `visualization::...`, `elliptic_curves::analytic::...`, or `numerics::...`.
- Treat `lib.rs` as a narrow convenience barrel, not as a mirror of the whole
  crate tree. Prefer re-exporting constructors, core traits, main strategy
  enums, and principal mathematical value types there, while keeping
  route-specific reports, step structs, and specialized diagnostics in their
  domain namespaces unless they are genuinely central to the crate's story.
- When deciding whether a symbol belongs in `lib.rs`, prefer the stricter
  question “would a first crate-level example naturally import this from the
  root?” over “could this save a longer path?”. If the main justification is
  path shortening for a specialized domain, keep it out of the crate root.
- Prefer namespace-first public APIs for this crate: `fields::Q`,
  `elliptic_curves::ShortWeierstrassCurve`, `polynomials::DensePolynomial`,
  and similar paths are the intended user-facing style. Do not reintroduce
  crate-root aliases for those values without a strong pedagogical reason.
- Prefer lightweight, mathematically honest type-level encodings when they
  remove the need for duplicate runtime context, as in `ExtensionField<S>`.
- Keep error ownership local to the domain:
  - `FieldError` in `fields`
  - `PolynomialError` in `polynomials`
  - `CurveError` in `elliptic_curves`
  - avoid duplicating the same failure mode as unrelated strings in several
    files
- Generic torsion-order logic belongs under `elliptic_curves`, not under
  `isogenies` or `division_polynomials`.
- Endomorphism-ring logic belongs under `elliptic_curves::endomorphisms`,
  not under `isogenies::graphs`. If graph helpers later consume that data,
  keep the arithmetic source of truth in `elliptic_curves` and treat any
  graph-side use as a consumer, not as the owner of ring logic.
- For the staged `PrimeNormIdeal` to binary-quadratic-form bridge, keep
  GP/PARI-derived fixtures under `elliptic_curves::endomorphisms` and treat
  `PrimeNormIdeal::root_mod_ell()` as a selected square root of the order
  discriminant modulo `ℓ`; do not promote the bridge to public API until tests
  fix whether the resulting form represents the ideal class or its inverse.
- In that crate-private bridge, choose the middle coefficient by the congruences
  `b ≡ root_mod_ell() (mod ℓ)` and `b ≡ Δ (mod 2)`, then build
  `(ℓ, b, (b² - Δ)/(4ℓ))` before Gauss reduction. Preserve both the raw and
  reduced forms in the internal report while later action-level APIs are still
  being designed.
- When any endomorphism-side algorithm has already chosen `a`, `b`, and a
  discriminant `Δ` for a binary quadratic form, use the shared crate-private
  `BinaryQuadraticForm::from_leading_middle_discriminant(...)` helper instead
  of recomputing `c = (b² - Δ)/(4a)` locally.
- When a caller only needs to know whether a reduced form belongs to one
  `QuadraticClassGroup`, prefer `contains_reduced_form(...)` over materializing
  and searching `enumerate_reduced_forms()` at the call site.
- When `quadratic_ideals` coverage grows, prefer a local `tests/` module tree
  split by testing intent, such as prime behavior, prime-norm ideal basics,
  construction errors, and GP/PARI ideal-form fixtures, instead of rebuilding a
  catch-all `tests.rs`.
- Division-polynomial-driven torsion search belongs under
  `elliptic_curves::short_weierstrass::division_polynomials`, even when later consumers are
  graph or isogeny features.
- For the first `endomorphisms` milestone, prefer a narrow finite-field
  surface that starts from existing Frobenius packages rather than claiming a
  full generic computation of `End(E)` or a certified maximal order.
  Good first surfaces are explicit reports or value objects derived from
  `FrobeniusTrace` / `FrobeniusCharacteristicPolynomial`; avoid introducing
  broader traits or graph-level classifications until the mathematics and
  tests already justify them.
- Division-polynomial explanation helpers and compact summaries belong under
  `visualization::elliptic_curves`, not under the core algebra modules.
- When a prime field is known at compile time, prefer the static Montgomery
  namespace `Fp<M, LIMBS>` with `ConstPrimeMontyParams`; keep `Fp<M, LIMBS>` as
  legacy small-field compatibility rather than the canonical direction for new
  APIs.
- When an algebraic extension can be described statically, prefer a
  specification type plus `ExtensionField<S>` so the extension still
  participates in the main `Field` trait and can itself serve as the base of a
  tower.
- When a higher tower step is mathematically valid but the crate does not yet
  have a generic irreducibility backend for that base field, a documented
  manual validation hook is acceptable in examples and educational extension
  specs. Mark that choice clearly as temporary.
- Do not smuggle ambient field context into element values when a cleaner
  field-family boundary is available.
- Avoid cross-module coupling unless it meaningfully improves clarity.
- Do not add new abstraction layers unless they remove real duplication or
  express a real mathematical boundary.
- When a capability is backend-specific, prefer a narrow trait such as
  `SqrtField`, `EnumerableFiniteField`, or a curve-side capability trait such
  as `LiftXCoordinate`, `EnumerableCurveModel`, `GroupCurveModel`, or
  `FiniteGroupCurveModel` over inflating a base trait that many backends
  cannot honestly implement.

## Development workflow

Before considering a change complete, run:

- `cargo fmt`
- `cargo test`
- `cargo clippy --all-targets --all-features`

Mutation-testing safety rule:

- Never stop halfway through a `cargo-mutants --in-place` rerun.
- If an in-place mutation run must be aborted because of an external failure,
  immediately inspect the touched file(s) for leftover
  `/* ~ changed by cargo-mutants ~ */` markers or equivalent mutated code and
  restore the original source before running anything else.
- Prefer non-`--in-place` mutation runs by default when the extra copy cost is
  acceptable; reserve `--in-place` for intentionally targeted reruns where the
  speedup is worth the extra operational risk.

When adding or modifying a runnable example, also run it once if it is cheap
and deterministic. The current dual-isogeny and graph examples should be
exercised with:

- `cargo run --example dual_isogeny`
- `cargo run --example isogeny_graph`
- `cargo run --example division_polynomials`

If a change is intentionally partial, the code should still compile and the
remaining work should be clearly signposted.

When a capability moves from “future work” to “implemented”,
update the top-level `README.md` in the same piece of work so the public repo
summary, example list, and feature narrative do not lag behind the code.
Keep that `README.md` high-signal and short: it should orient a reader to the
project, not try to enumerate every implemented surface or become a changelog.

When changing the control flow, case split, pipeline stages, or key invariants
of an algorithm that is documented with Mermaid diagrams under
`docs/algorithm-diagrams/`, update the corresponding diagram in the same piece
of work. Treat those diagrams as educational API documentation: stale
diagrams are a correctness and pedagogy bug, not just a docs nit.

### Verification strategy for new features

When advancing the repository with new mathematical features, prefer this
verification ladder:

1. start with Rust tests
2. add small exhaustive searches in Rust or a sidecar when the state space is
   still honestly tiny
3. use SMT, with cvc5 as the default solver, as an auditor of discrete
   invariants once the feature makes a strong algebraic claim or crosses a
   mathematically subtle boundary
4. reserve heavier formalization work such as Lean for central, stable
   mathematical statements that have already proved their value in the codebase

In particular, when a feature introduces a claim like:

- “same invariant implies same structure”
- “every candidate produced by this algebraic filter is exact”
- “this exceptional family behaves just like the generic one”

prefer asking cvc5 to search for a small counterexample before treating the
claim as trusted. Use SMT primarily as a bug hunter first and only secondarily
as a proof-style checker once the search for violating models has failed.

## Testing rules

- Do not add large algorithms without tests.
- Prefer deterministic, small examples first.
- For algebraic structures, add property-oriented tests where appropriate:
  - associativity
  - identity laws
  - inverses
  - distributivity
  - compatibility with canonical reduction
- For capability traits such as square roots or curve membership, test both the
  success path and the honest “no solution / not supported” path that the API
  promises.
- For educational helpers such as formatting or visualization functions, test
  the textual output at the level of important content, not brittle full-file
  snapshots unless the output format is intentionally fixed.
- When a module exposes typed errors, test the error variants directly instead
  of asserting only on formatted messages.
- When polynomial or quotient representations are added, include tests for both
  data invariants and how the chosen storage order is explained to readers.

## Performance rules

- Do not optimize early.
- Prefer the clearest correct implementation first.
- If performance-oriented code is added later, preserve a readable reference
  path when possible.
- Avoid introducing specialized arithmetic tricks, unsafe code, or heavy
  dependencies without a concrete demonstrated need.

## Dependency policy

- No dependency should be added casually.
- A new dependency must have a narrow, justified purpose.
- If a dependency is added, keep the usage small and document why it belongs.
- Prefer standard library facilities unless an external crate materially
  improves correctness, clarity, or educational value.

Current justified numeric dependencies include:

- `num-complex` for approximate complex arithmetic
- `num-bigint` and `num-rational` for exact arithmetic over `Q`
- `num-traits` for numeric identities used by those exact types

For small exact integer infrastructure already backed by `num-bigint`, prefer
one canonical helper in `numerics` for reusable operations such as `gcd` and
`lcm` on `BigUint` instead of scattering local copies across domains.

## Scope discipline

- This repository is not yet a production cryptography crate.
- Do not present scaffold code as production-safe.
- Do not harden APIs prematurely around features that are not implemented yet.
- Avoid speculative support for serialization, randomness, parallelism, or FFI
  unless the project explicitly moves in that direction.
- Do not assume every algebraic construction should be phrased as a finite
  field. Infinite fields such as `Q` are first-class educational citizens in
  this codebase.

## Final reporting expectations

When summarizing work:

- mention the main files changed
- describe the conceptual change, not just the diff
- note any important simplifications made
- mention remaining risks, TODOs, or intentionally deferred work

The best changes in this repository should feel mathematically honest,
pedagogically useful, and easy for the next contributor to continue.
