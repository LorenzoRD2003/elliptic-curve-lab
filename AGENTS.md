# AGENTS.md

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
- Support both finite and infinite base fields when the mathematics naturally
  calls for it, instead of assuming everything is cryptographic or finite from
  the start.
- Expose mathematically meaningful backend properties when they improve later
  APIs, such as whether a field family is algebraically closed.
- Keep capability boundaries explicit when only some backends honestly support
  an operation, as with square roots or future curve-side helpers.
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
    `elliptic_curves::division_polynomials`
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

At the moment, the most mature parts of the repository are `fields` and
`polynomials`, especially:

- `Fp<P>` and `FpElem<P>` for exact prime-field arithmetic
- `Q` for exact rational arithmetic over `BigRational`
- `ComplexApprox` for approximate numerical experiments over `C`
- `SqrtField` as a small capability trait for backends that can produce square
  roots honestly, with current implementations for `Fp<P>`, `Q`, and
  `ComplexApprox`
  - a quadratic-character capability for finite fields of odd characteristic,
  with explicit `0 / residue / non-residue` values rather than a hidden
  ad hoc integer convention
  - a quadratic-character point-count route
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
- the first usable pieces of `elliptic_curves`, currently centered on affine
  points, short-Weierstrass curves, discriminants, curve-membership checks,
  `x`-coordinate lifting, small-field point enumeration, a first explicit
  group-law trait for additive curve operations, small-group helpers such as
  torsion checks and point orders, classical short-Weierstrass invariants
  such as `c4`, `c6`, and `j`, plus a first explicit Frobenius layer that
  distinguishes:
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
  - the next planned algebraic step is an educational
    `elliptic_curves::endomorphisms` layer derived from finite-field
    endomorphism-ring information already visible through Frobenius-side data
    such as `q`, `t`, `χ_{π_q}(T)`, and `t^2 - 4q`
  - a first function-field layer for short-Weierstrass curves, modeling
    `F(E) = F(x) ⊕ yF(x)` through pairs of rational functions
    `(A(x), B(x))` representing `A(x) + yB(x)`, with multiplication reduced by
    the specific short-Weierstrass relation `y^2 = x^3 + ax + b`, plus public
    substitution helpers for evaluating polynomials and rational functions in
    the distinguished `x`-coordinate at a function-field element
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
- the first usable pieces of `isogenies`, including explicit finite kernels,
  Vélu isogenies on short-Weierstrass curves, exhaustive structural
  verification helpers, strict and bridged composition on small finite curves,
  scalar-multiplication isogenies `[n]`, exhaustive map-equality helpers, and
  dual Vélu search by enumerating tiny rational kernels and testing both
  duality relations on rational points, plus public helpers for checking
  `\hat{\phi} \circ \phi = [deg(\phi)]` and
  `\phi \circ \hat{\phi} = [deg(\phi)]` exhaustively, together with a first
  short-Weierstrass pullback layer on function fields that represents
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
- text-based visualization helpers for dual-isogeny workflows,
  including composition summaries, scalar-multiplication summaries, dual
  isogeny summaries, and exhaustive dual-verification reports suitable for the
  final dual-isogeny example
- runnable educational examples under `examples/`, including extension towers
  plus walkthroughs for curve order, group structure, isomorphisms,
  Vélu isogenies, dual isogenies, `ℓ`-isogeny-graph exploration, and
  division-polynomial torsion recovery that show how the APIs and
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
  math. Do not use backticks such as ``...`` or fenced blocks like
  ```math``` for mathematical notation.
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
- Division-polynomial-driven torsion search belongs under
  `elliptic_curves::division_polynomials`, even when later consumers are
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
- When a field family is known at compile time, prefer a namespace type such as
  `Fp<P>`.
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
