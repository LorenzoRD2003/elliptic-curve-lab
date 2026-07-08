use crate::fields::traits::*;
# AGENTS.md for `src/elliptic_curves`

## Module mission

The `elliptic_curves` module should introduce curve models gradually and
honestly.

Right now the goal is not to ship a production EC library. The goal is to make
curve structure, point representations, and equation checks easy to read and
easy to extend.

## Current posture

- Early-stage scaffolding is acceptable when it is explicit and tested.
- Short Weierstrass support is currently the main concrete path.
- For staged work on new curve families such as `general_weierstrass`, keep the
  implementation slices narrow and add local tests as each slice lands instead
  of deferring all coverage to the final integration pass.
- During local verification for edits under `src/elliptic_curves`, prefer
  `cargo test -q` with a module-level filter that matches the code being
  changed, then broaden only when the edited surface actually crosses module
  boundaries. Avoid the full `cargo test -q` for localized work in this tree
  unless the scope genuinely spans enough neighboring modules that targeted
  verification would stop being honest.
- At the start of each new conversation affecting this module tree, first
  consult an existing Repomix artifact if one is present, or generate a
  compressed one if not, and use its token summary to decide whether large
  guidance files, long docs, or oversized source files should be excluded from
  immediate context before deeper analysis begins.
- For localized curve-model work, keep immediate context narrow by default:
  the repo-root `AGENTS.md`, this module-local `AGENTS.md`, the active model
  plan, and nearby code should stay in focus first. Avoid loading unrelated
  heavyweight files such as distant visualization/Frobenius surfaces or broad
  algorithm-diagram docs unless the active task genuinely depends on them.
- For reductions from one curve model to another, prefer an explicit reduction
  value object that stores the source model, the target companion, and the
  coordinate-change parameters, rather than exposing only the target curve.
- If one reduction direction is just the obvious embedding, still prefer
  encoding it explicitly in that same reduction layer so roundtrip tests and
  later trait integrations have one symmetric place to depend on.
- For public cross-model APIs, prefer one reusable conversion trait in
  `elliptic_curves::traits` and keep concrete reduction witnesses private to
  the owning model family unless callers genuinely need the raw parameters.
- Do not reintroduce an `ambient_short_weierstrass` lane or a `BigPrimeField`
  runtime curve stack. That experiment duplicated the existing curve logic too
  much. `AmbientField` may still appear inside function-field implementation
  bodies, but it should not drive a public curve family or a Cargo feature.
- For the rational-torsion reduction-mod-`p` route, keep runtime-prime
  arithmetic local to `short_weierstrass::rational_torsion::reduction_mod_p`
  with role-specific names such as `ReductionPrime`/`ReductionResidue`; do not
  promote it into `fields` or make it implement the public `Field` traits.
- In that same rational-torsion reduction route, the first educational
  good-prime search should start at `p = 11`, skip primes dividing the integral
  discriminant `Δ`, and record the nonzero residue `Δ mod p` as the local
  certificate of good reduction before later point-enumeration stages use it.
- For the first reduced-curve stage of that route, keep `E_p(𝔽_p)` enumeration
  deliberately exhaustive and local to `reduction_mod_p`: reduce `A,B` through
  the runtime-prime helper, scan affine pairs in `Θ(p²)`, and only optimize
  after a later stage has a measured need or broader reuse.
- For the reduced-curve Mazur filter, use the same fixed point-order list as
  the rational exact verifier and keep the result as reduced candidates only:
  filtering `E_p(𝔽_p)` by `[m]P = O` for Mazur-permitted `m` is not yet a
  rational lift or a torsion classification until the later Hensel/exact
  verification stages certify it.
- When the rational-torsion reduced curve needs short-Weierstrass group-law
  formulas, adapt `group_law_core` through local runtime-residue operations
  rather than duplicating slope/reconstruction formulas inside
  `reduction_mod_p`.
- When the rational-torsion reduction route needs order-specific
  `x`-criteria, reuse `models::short_weierstrass::division_polynomials` as
  the canonical source for `ψ_m` and the even-index `y`-factor split. The
  `reduction_mod_p` layer may adapt those criteria to primitive
  `IntegerPolynomial` values for modular seeds, but should not rederive or
  duplicate division-polynomial recurrences.
- Keep that `x`-criterion adapter scoped to the Mazur-permitted non-identity
  point orders used by rational torsion over `Q`. Do not let it silently
  become a general-purpose division-polynomial wrapper for arbitrary indices.
- In that same route, for even `m` it is acceptable to use the stored factor
  `f_m(x)` from `ψ_m = y f_m(x)` instead of literally forming `ψ_m/ψ_2`,
  because `f_m = 2(ψ_m/ψ_2)` and the scalar `2` is invertible over `Q` and
  modulo the chosen primes `p ≥ 11`. Represent that scalar-equivalence with a
  distinct adapter source variant, and document the primary mathematical
  criterion as `ψ_m/ψ_2` rather than as the implementation factor `f_m`.
- For the first Hensel stage of the rational-torsion reduction route, lift
  only the `x`-coordinate seeds produced by `TorsionXPolynomial`, use the
  polynomial-side Cauchy bound for integer-root recovery, and record singular
  or uncertified seeds as report outcomes. Do not jump directly from modular
  seeds to rational points until a later stage performs `y`-recovery and exact
  torsion verification.
- Once that route performs `y`-recovery, keep it on the integral companion
  model first: solve `y² = x³ + Ax + B` exactly over `ℤ`, construct rational
  points on the integral curve, and verify them by exact Mazur-order scalar
  multiplication. Leave transport back to the source curve and Lutz-Nagell
  comparison to a later integration/reporting pass.
- Keep the first point-lift report small: store the chosen good prime, the
  verified group, and the verified points. Do not duplicate polynomial/Hensel
  subreports or discard counters in that layer until the educational report
  format explicitly needs them.
- For public rational-torsion computation over `Q`, keep
  `ShortWeierstrassCurve::rational_torsion_by(RationalTorsionStrategy)` as the
  single caller-facing entry point. Route-specific engines such as
  Lutz-Nagell enumeration and good-reduction/Hensel lifting should remain
  internal and feed the shared `RationalTorsionReport`.
- Keep rational-torsion timing comparisons as ignored tests or external
  benchmarks, with a deterministic corpus and an explicit correctness check
  before timing. Normal module tests may compare strategy outputs, but should
  not assert wall-clock speed.
- When a helper in the point-lift route only sorts/deduplicates verified
  `(point, order)` pairs, keep its input and output representation the same.
  Perform `unzip` or other shape changes at the call site where the next
  consumer actually needs separate vectors. After sorting by point, prefer
  `dedup_by` over a hand-rolled membership loop.
- In Hensel-stage reports, prefer storing the successful lift trace as the
  canonical payload and deriving certified integer roots from that trace rather
  than caching a duplicate root field alongside it.
- For complexity rustdocs inside Mazur-bounded rational-torsion helpers, treat
  the fixed list of permitted orders as constant-sized and state costs in
  terms of visible variable input sizes such as coefficient count or prime
  size. Prefer simple expressions like `Θ(p·n)` over opaque placeholders for
  nested subroutine costs when the precise expansion would obscure the
  algorithmic story.
- Do not treat `F::characteristic()`, scalar multiplication, Frobenius traces,
  Hasse intervals, or isogeny degrees as fixed-width integers in new curve
  code. Use `F::has_characteristic(2)` / `F::has_characteristic(3)` for model
  restrictions, and otherwise carry `BigUint`/`BigInt` through curve reports
  and public educational APIs. If an executable search engine still needs a
  primitive index because it materializes a table or loops over an enumerable
  tiny set, keep that primitive local to the engine and do not let it leak into
  report structs or scalar APIs.
- For the first quadratic-ideal layer, keep
  `endomorphisms::quadratic_ideals` focused on local prime behavior in an
  imaginary quadratic order. The public `prime_behavior(ℓ)` surface may answer
  split/inert/ramified/non-invertible behavior for `(Δ/ℓ)`, but it should not
  introduce ideal objects, ideal classes, composition, or curve actions until
  those later layers have their own documented invariants.
- When that layer introduces prime-norm ideal objects, keep the first public
  surface to an opaque `PrimeNormIdeal` with crate-internal split-prime
  and ramified-prime representations: local root modulo `ℓ`, norm,
  conjugation, and validation are enough. Do not add a shared trait, ideal
  multiplication, class equivalence, or curve-side kernels/actions in the same
  step.
- Examples for complex analytic curves should require the `analytic` Cargo
  feature, while examples for Schoof, Mestre, or Hasse-search comparison
  routes should require `advanced-point-counting`. These feature names mark
  educational chapters first; do not split them into per-algorithm flags until
  there is a concrete compile-time payoff.
- When a curve-side helper is used only by visualization tests, gate it with
  `#[cfg(all(test, feature = "visualization"))]` so `--all-features --examples`
  does not compile test-only helpers as dead code.
- If such a shared trait requires one small foundational capability from a
  staged model, such as `CurveModel`, prefer implementing just that minimal
  capability now instead of weakening the trait contract globally.
- When a model-side trait introduces its own error surface, prefer wiring the
  shared and local error types together with `From` where the conversion is
  total, and reserve handwritten `map_err` matches for genuinely contextual
  cases such as source-vs-target point failures.
- If one internal helper error for a concrete curve family always collapses to
  `CurveError` without extra context, prefer `From<HelperError> for CurveError`
  over a local mapper function living next to one trait impl.
- When an explicit conversion witness already exists for two curve families,
  prefer also exposing whole-curve `From` or `TryFrom` impls for the
  companion-model conversion itself, so callers can reuse the same model
  relationship without rebuilding ad hoc helper names.
- Treat `LiftXCoordinate` as the finite affine-fiber story above one chosen
  `x`, not as a short-Weierstrass-only `y^2 = rhs(x)` helper. Different curve
  families may realize that fiber through square roots, shifted quadratics, or
  future Artin-Schreier solvers, and the trait should stay honest about that.
- As a stepping stone toward general-Weierstrass lifting, prefer factoring the
  `y^2 + uy = v` equation into its own small helper with explicit odd-characteristic
  and future characteristic-`2` paths, rather than baking that algebra into one
  premature public trait.
- Once the field layer exposes characteristic-`2` Artin-Schreier solving,
  prefer wiring that helper all the way through `GeneralWeierstrassCurve`'s
  `LiftXCoordinate` implementation so the trait talks about actual `x`-fibers
  instead of leaving the characteristic-`2` branch as a documented gap.
- Within that general-Weierstrass fiber helper, prefer names that describe the
  geometry directly, such as `y_fiber_equation`, `linear_coefficient`, and
  `right_hand_side`, and keep `model_traits.rs` as a thin adapter from solved
  `y`-fibers into checked affine points.
- If that same helper grows enough to mix equation state, backend solvers, and
  curve-side extension methods, prefer a `y_fiber/` submodule with focused
  siblings such as `equation.rs`, `solvers.rs`, and `curve.rs`.
- Once that split exists, prefer moving tests whose subject is really the
  `y_fiber` helper itself into `y_fiber/tests.rs`, while leaving broader
  curve-model or trait-integration tests at the parent `general_weierstrass`
  level.
- When `GeneralWeierstrassCurve<F>` gains `GroupCurveModel` support in staged
  form, prefer starting with honest native affine formulas, including the
  model-specific negation involution, before moving on to projective formulas.
- For staged Montgomery work, do not reject characteristic `3` at model
  construction time: the classical affine Montgomery model only requires
  characteristic different from `2`, while the familiar short-Weierstrass
  companion conversion is the step that needs separate characteristic-`3`
  handling.
- When implementing Montgomery invariants, cross-check them either against the
  equivalent general-Weierstrass model or against the identities
  `j = c4^3 / Δ` and `c4^3 - c6^2 = 1728Δ`; the intended normalization here is
  `c4 = 16(A^2 - 3)/B^2`, `c6 = 32A(9 - 2A^2)/B^3`, and
  `Δ = 16(A^2 - 4)/B^6`.
- For staged Montgomery `LiftXCoordinate`, prefer the direct fiber equation
  `y^2 = (x^3 + A x^2 + x)/B` with a `SqrtField` bound; once that is in place,
  let `EnumerableCurveModel` come for free from the blanket trait rather than
  building a model-specific enumeration surface first.
- For staged Montgomery reduction to `ShortWeierstrassCurve<F>`, prefer one
  explicit witness that stores the Montgomery source, the short companion, and
  the affine transport determined by `x = B X - A/3`, `y = B Y`; keep this
  route unavailable in characteristic `3` even though the Montgomery model
  itself remains valid there.
- For staged Montgomery-to-general conversion, prefer the direct whole-curve
  embedding with coefficients `a1 = 0`, `a2 = A/B`, `a3 = 0`, `a4 = 1/B^2`,
  `a6 = 0` instead of routing through the short companion when the task only
  needs a general-Weierstrass view of the curve.
- For staged short/general-to-Montgomery conversion without a user-supplied
  witness, it is acceptable to provide only a curve-level `TryFrom` route
  under finite enumerable plus square-root-capable bounds, driven by a
  rational `2`-torsion root and the tangent factor `3α^2 + a`; keep failure
  explicit when that certification does not exist over the current base field.
- For staged Montgomery `GroupCurveModel` support, prefer one honest native
  affine implementation before any projective layer: negation is
  `(x,y) -> (x,-y)`, secant addition for distinct `x` uses
  `λ = (y2-y1)/(x2-x1)` with
  `x3 = B λ^2 - A - x1 - x2`,
  `y3 = λ(x1-x3) - y1`, and doubling uses
  `λ = (3x^2 + 2Ax + 1)/(2By)` with
  `x([2]P) = B λ^2 - A - 2x`,
  `y([2]P) = λ(x-x([2]P)) - y`.
- Once that staged affine Montgomery group law exists, validate it against the
  short companion in characteristic `> 3` and with exhaustive small-field
  group-axiom checks in at least characteristics `3` and `5` before moving on.
- For staged Montgomery ladder work, treat the `x`-coordinate story as its own
  milestone: introduce explicit Montgomery-owned `x`/`X:Z` value objects first,
  then documented differential primitives such as `xDBL` and `xADD`, and only
  then the scalar ladder itself. Do not hide the first ladder implementation as
  “constant-time” unless that stronger claim is separately justified for the
  actual backend and control-flow details.
- For the first Montgomery ladder implementation in this repo, prefer an
  explicit internal normalization to the classical `B = 1` ladder formulas
  when, and only when, the required scaling witness exists over the current
  base field; if that witness is unavailable, fail honestly rather than
  pretending the normalized ladder applies to every `B y^2 = x^3 + A x^2 + x`
  model uniformly.
- The intended first normalization witness for Montgomery ladder work is the
  `y`-rescaling `v = sqrt(B) y`, so normalization over the same base field
  exists exactly when `B` is a square there. Prefer recording that witness as a
  first-class reduction-like object that stores the source curve, the
  normalized `B = 1` target, and the chosen `sqrt(B)`, and treat the
  non-square case honestly as “no same-field normalization witness” rather than
  as an invisible implementation detail.
- For the staged Montgomery `x`-line representation, prefer an explicit
  infinity variant plus finite projective `X:Z` representatives up to scaling;
  affine `x = X/Z` recovery and normalization should require `Z` to be
  invertible rather than silently overloading `Z = 0` as a second infinity
  encoding.
- For the first normalized Montgomery differential-arithmetic layer on
  `v^2 = x^3 + A x^2 + x`, use the classical `X:Z` formulas with
  `A24 = (A + 2)/4`: doubling should follow
  `AA = (X+Z)^2`, `BB = (X-Z)^2`, `E = AA-BB`,
  `X([2]P) = AA*BB`, `Z([2]P) = E*(BB + A24*E)`, and differential addition
  should follow
  `DA = (XP-ZP)(XQ+ZQ)`, `CB = (XP+ZP)(XQ-ZQ)`,
  `X(P+Q) = Z(P-Q)(DA+CB)^2`, `Z(P+Q) = X(P-Q)(DA-CB)^2`.
- For the first Montgomery ladder stage, keep the maintained state as one
  neighboring pair `(R0, R1)` with fixed difference `R1 - R0 = P`, initialize
  it as `(O, P)`, and update it only through the Stage B `xDBLADD` pattern.
  Document that this is a fixed-schedule educational ladder, not yet a
  production constant-time claim across all field backends.
- When documenting the current Montgomery ladder layer, state its computational
  complexity explicitly: `Θ(log n)` differential steps for an arbitrary-size
  scalar `n`, with one uniform `xDBLADD`-shaped schedule step per processed
  bit, and keep that complexity note in rustdocs near the executable ladder
  entry points.
- Once the Montgomery ladder exists, add both exhaustive tiny-field validation
  in representative characteristics such as `3` and `5` and broader property
  tests on normalized cases; keep the public explanation explicit that the
  ladder returns an `x`-coordinate class, not a signed affine point, and prefer
  one small report type when that makes the limitation visible to callers.
- For staged Montgomery finite-field APIs, keep one curve-side wrapper per
  invariant family, matching the existing general-Weierstrass story:
  `group_order_by(...)`, `group_order_by_small_field(...)`,
  `point_order_by(...)`, and `group_exponent_by(...)`. Let `Exhaustive` and
  other purely enumerative routes stay native in characteristic `3`, but keep
  quadratic-character and Schoof-style routes delegated through the short
  companion only when the reduction to short-Weierstrass is available.
- For early CM trace-candidate work, keep the arithmetic helper under
  `frobenius::cm`: it may use quadratic forms and Cornacchia internally, but
  it should report only Frobenius-side candidates such as `|t|` until a later
  curve-specific API certifies that a concrete curve has CM by the supplied
  discriminant and determines the sign of the trace.
- For staged `TwistedEdwardsCurve<F>` work, prefer characteristic different
  from `2` as the first milestone, with honest descriptor validation through
  `a != 0`, `d != 0`, and `a != d` before any deeper executable layer lands.
- For that same Twisted-Edwards Stage A descriptor milestone, prefer keeping
  the executable surface minimal: validated coefficients, equation formatting,
  and classical invariants. If the invariant formulas are derived through the
  canonical Montgomery bridge normalization, document that normalization
  explicitly near the implementation instead of hiding it behind unexplained
  constants.
- For that same Twisted-Edwards family, prefer making `MontgomeryCurve<F>` the
  canonical whole-curve bridge: own the direct coefficient formulas at the
  Edwards/Montgomery boundary, then reuse existing Montgomery-to-short and
  Montgomery-to-general routes by composition instead of adding parallel
  Twisted-Edwards reductions immediately.
- Treat Twisted-Edwards point transport to Montgomery as a separate milestone
  from whole-curve conversion. The standard formulas are birational and have
  exceptional affine loci, so do not promise the current total-point
  `CurveModelConversion` contract on day one unless the chosen point
  representation or bridge abstraction really closes those exceptional cases.
- Twisted Edwards is also the first likely affine family here with a finite
  neutral element `(0, 1)`. Before relying on blanket enumerable helpers,
  prefer making the shared enumeration path identity-aware rather than
  encoding the neutral element artificially as the point at infinity.
- In that shared enumerable path, `finite_points()` should filter out the
  model-defined identity even when `lift_x` returns it as a finite affine
  point, and `points()` should then add the identity back exactly once at the
  front of the returned list.
- For staged Twisted-Edwards transport APIs, prefer explicit birational-open
  naming such as `try_point_to_montgomery_open(...)` over names that suggest a
  total affine point equivalence.
- For the Stage C whole-curve Twisted-Edwards/Montgomery bridge, keep the API
  descriptor-only and total: coefficient-level conversion is honest here, so
  prefer `as_montgomery()` / `as_twisted_edwards()` and `From<&...>` impls,
  while continuing to reserve point-transport APIs for the later birational
  stage.
- Keep the semantic distinction explicit between:
  - whole-curve conversion
  - birational point transport on an affine open subset
  - later total rational-point correspondence
  especially before adding ladders, torsion/cofactor surfaces, or isogeny-side
  transport that could otherwise inherit the wrong mental model.
- For the first native affine `TwistedEdwardsCurve<F>` group law, prefer the
  generic `(a, d)` family with honest denominator checks and documented affine
  formulas; reserve "complete formula" claims for a later, separately scoped
  restricted-subfamily milestone.
- For the first Twisted-Edwards projective step, prefer a local
  `twisted_edwards/projective/` directory and make the extended-point layer
  explicit before any native projective group law lands: use
  `ExtendedTwistedEdwardsPoint<F>` in `(X:Y:Z:T)` coordinates, validate it by
  `aX^2 + Y^2 = Z^2 + dT^2` together with `XY = ZT`, keep the affine embedding
  `(x, y) -> (x : y : 1 : xy)`, and let `to_affine` fail honestly when
  `Z = 0` instead of inventing a hidden affine fallback.
- For Twisted-Edwards test layout, keep implementation helpers under
  `twisted_edwards/projective/` but place their coverage in the family-level
  `twisted_edwards/tests/` tree, grouped by responsibility, rather than
  introducing a second nested test hierarchy under the implementation module.
- When Twisted Edwards gains `HasProjectiveModel` before a native projective
  group law, wire that trait directly to `ExtendedTwistedEdwardsPoint<F>` and
  keep the projective identity equal to the finite neutral representative
  `(0:1:1:0)`, not to a Weierstrass-style infinity sentinel.
- For the first `ProjectiveGroupCurveModel` milestone on Twisted Edwards,
  prefer an explicit affine-bridge implementation over premature native
  extended-coordinate formulas: validate extended-point membership, recover to
  affine, delegate to the existing affine group law, and lift back. Keep that
  bridge documented as transitional until the extended formulas are
  independently implemented and validated.
- Once native extended-coordinate addition and doubling land for Twisted
  Edwards, prefer promoting that projective engine to the single executable
  group-law core: let the affine public wrappers delegate through
  `to_projective -> native projective op -> to_affine_projective`, and keep
  any resulting `Z = 0` affine-recovery failure as the honest way generic
  affine denominator singularities surface.
- When documenting the twisted-Edwards extended-coordinate group law, state
  explicitly that `Z = 0` means the result still belongs to the projective
  model but has left the affine chart `x = X/Z`, `y = Y/Z`; the public affine
  wrappers should then be described as failing at the recovery step, not as if
  the projective group law itself had become undefined.
- In that first staged affine Twisted-Edwards group law, when a valid
  on-curve input hits a zero denominator in the generic formulas, prefer
  surfacing that honestly as `CurveError::Field(FieldError::DivisionByZero)`
  rather than degrading it to `PointNotOnCurve` or silently branching into a
  stronger completeness claim than the current formulas actually justify.
- For the Stage F finite-compatibility milestone on `TwistedEdwardsCurve<F>`,
  prefer relying first on the shared `EnumerableCurveModel` and
  `FiniteGroupCurveModel` surfaces rather than adding new model-owned wrapper
  methods prematurely. Close the milestone with explicit tests that the finite
  identity is enumerated exactly once and that order/exponent/group-structure
  data agree with the Montgomery companion when comparison is mathematically
  honest.
- For staged Twisted-Edwards membership and `CurveModel` support, treat the
  neutral element `(0, 1)` as the only identity point and reject
  `AffinePoint::Infinity` as not belonging to the affine model. The identity
  is finite here because that is part of the model's actual geometry, not a
  transport artifact.
- When Montgomery reaches the educational/examples milestone, prefer one
  runnable example and one visualization helper set that show the native
  Montgomery equation, the short companion when available, the direct general
  embedding, and at least one point-level or group-law comparison transported
  across those models.
- Keep an explicit TODO next to that affine implementation that the intended
  long-term replacement is a projective-coordinate general-Weierstrass group
  law.
- When validating that staged group law, prefer exhaustive checks over tiny
  finite fields in characteristics such as `2` and `3`, and back them with
  reusable `proptest_support::elliptic_curves` strategies for broader
  property coverage.
- When integrating the static `crypto-bigint` prime-field backend with curve
  families, test construction, membership, `LiftXCoordinate`, and existing
  group/projective scalar surfaces directly over `Fp<M, LIMBS>` without
  adding `EnumerableFiniteField` bounds or reusing exhaustive small-field
  helpers.
- Once a new curve family reaches blanket compatibility with finite-group or
  Frobenius-side traits, add dedicated compatibility tests for enumeration,
  point orders, group structure, Frobenius trace/Hasse workflows, and cyclic
  kernel construction instead of assuming the blanket integration is enough by
  inspection alone.
- Prefer organizing those tests under a local `tests/` directory split by
  intention, following the existing short-Weierstrass style, rather than
  keeping one mixed `tests.rs` once the family has several independent test
  stories.
- For staged ergonomic wrappers on `GeneralWeierstrassCurve<F>`, prefer
  exposing the familiar curve-side APIs only when their docs say explicitly
  which routes are native, which routes are delegated to the short companion,
  and which characteristics remain outside the delegated routes.
- For `GeneralWeierstrassCurve<F>`, keep transport-certification tests in a
  dedicated `tests/compatibility.rs` file once the model coexists with the
  short companion, so invariants, point transport, and order/group-law
  compatibility stay easy to audit as a single contract.
- Once that same `GeneralWeierstrassCurve<F>` milestone reaches the
  “deep generalization” planning stage, keep the long-horizon plan as a
  versioned `roadmap.md` inside the model directory itself, and organize it by
  dependency order plus exit criteria rather than as one flat wishlist.
- For the odd-characteristic side of that same `GeneralWeierstrassCurve`
  lifting story, prefer relaxing the bounds back down to `Field + SqrtField`
  so exact characteristic-`0` backends such as `Q` and approximate complex
  backends such as `ComplexApprox` can participate before any finite-field-only
  machinery is involved.
- If rational extension fields are still excluded from that route, say so
  explicitly in docs: the blocker is the absence of a generic square-root
  backend there, not any obstruction in the `x`-fiber formulation itself.
- The top-level `elliptic_curves` tree should now stay organized around
  explicit ownership boundaries:
  - `affine.rs` for the shared affine point representation
  - `error.rs` for shared curve-domain errors
  - `models/` for concrete curve families, capability traits, and
    short-Weierstrass-specific executable layers such as division polynomials,
    function fields, isogenies, and isomorphisms
  - `group_algorithms/` for shared additive-group/order scaffolding that is not
    itself owned by one model family's mathematical definition
  - specialized theories such as `frobenius/`, `analytic/`, and
    `endomorphisms/` should remain their own namespaces instead of being
    flattened back into one top-level grab bag
- The public root of `elliptic_curves` should stay intentionally austere:
  - root items should be limited to core anchors such as `CurveError`,
    `AffinePoint`, and `ShortWeierstrassCurve`
  - model traits should be reached through `elliptic_curves::traits::<Trait>`
  - specialized reports, strategy enums, analytic values, Frobenius helpers,
    endomorphism-side types, function-field values, and isomorphism helpers
    should be reached through their own namespaces, not re-exported from the
    broadest barrel
- The first function-field layer should stay short-Weierstrass-specific under
  `models::short_weierstrass::function_fields` until the repo has a genuinely
  justified abstraction for other curve models.
- Public substitution helpers that are intrinsic to the short-Weierstrass
  function field, such as evaluating polynomials or rational functions in the
  distinguished `x`-coordinate at a function-field element, belong under
  `elliptic_curves::short_weierstrass::function_fields` rather than under
  downstream consumers such as `isogenies`.
- Likewise, generic-point arithmetic in `E(F(E))` belongs under
  `elliptic_curves::short_weierstrass::function_fields`: if Vélu, scalar
  multiplication, or other consumers need rational addition, translation, or
  doubling formulas on the generic point, prefer exposing one shared point
  representation and helpers there instead of re-encoding secant/tangent
  formulas in each consumer.
- The current affine representation should preserve mathematical invariants in
  the type when possible.
- Validation logic such as discriminant checks and point-membership checks is
  part of the educational API surface, not incidental glue.
- Classical short-Weierstrass invariants such as `Δ`, `c4`, `c6`, and
  `j` are appropriate here when their docs explain the mathematics directly.
- Curve-side capability traits are now part of the intended architecture:
  - `AffineCurveModel` for checked affine construction
  - `GroupCurveModel` for models that expose actual point addition, doubling,
    and scalar multiplication
  - `LiftXCoordinate` for models that can recover points from `x`
  - `EnumerableCurveModel` only for small finite settings where exhaustive
    point listing is honest
  - `FiniteGroupCurveModel` only for small finite settings where point orders
    and related order-theoretic helpers can be computed by direct group traversal
  - invariant capability traits such as `HasJInvariant` when a family can
    expose classical invariants honestly without inflating `CurveModel`
  - those invariant capability traits belong under
    `elliptic_curves::traits`, not under a separate `core` or
    `invariants` namespace
- A small `CurveIsomorphism` trait is now part of the intended architecture for
  explicit curve-to-curve witnesses. It should stay narrow: domain, codomain,
  and point evaluation.
- Exhaustive base-field isomorphism search on enumerable finite fields is now
  a first-class educational tool, not just an internal convenience. It is
  acceptable to use it to support dual-isogeny workflows such as dual-isogeny
  search, provided the docs say clearly that this is a tiny-field exhaustive
  routine.
- Division-polynomial and torsion helpers are now part of the
  intended `elliptic_curves` surface. It is acceptable to:
  - expose low-degree explicit division polynomials first
  - use recursive formulas plus memoization over small fields
  - compare division-polynomial torsion recovery against exhaustive point
    enumeration when the docs say so directly
  - prefer making those division-polynomial and torsion routines methods on
    `ShortWeierstrassCurve<F>` rather than public free functions in the module
  - when that layer grows, prefer splitting it by responsibility such as
    types, construction, evaluation, torsion search, and tests, instead of
    keeping one large mixed file
  - for point-order recovery from a known annihilating multiple, prefer a
    dedicated short-Weierstrass helper and a report that shows the per-prime
    peeling steps, rather than collapsing the algorithm to a bare integer
  - if that same algorithm factors into one reusable local `ℓ`-primary group
    routine plus one model-specific wrapper that isolates each prime
    component, prefer making the local additive-group routine the internal
    source of truth
  - but if that same workflow also needs exact integer infrastructure such as
    normalized prime-power factorizations or cached powers `1, ℓ, ..., ℓ^e`,
    prefer keeping those helpers in `numerics` and importing them here rather
    than rebuilding the arithmetic locally
  - the implementation owner of short-Weierstrass division-polynomial logic
    is now `models::short_weierstrass::division_polynomials`; do not recreate
    a sibling top-level `elliptic_curves::division_polynomials` tree
- Frobenius helpers are now part of the intended `elliptic_curves` surface.
  It is acceptable to:
  - expose the absolute Frobenius `π_p` separately from the relative
    Frobenius `π_q`
  - expose Frobenius-trace data as the trace of the relative Frobenius `π_q`
    of the chosen finite base field `F_q`, even when the current computation
    is still based on exhaustive point counting
  - prefer the notation `π_p` / `π_q` in mathematical docs, examples, and
    explanatory text instead of `Frob_p` / `Frob_q`
  - prefer associated constructors such as
    `AbsoluteFrobenius::for_field::<F>(...)` and
    `RelativeFrobenius::for_field::<F>(...)` over free module-level builder
    functions when the value being built is metadata for a concrete type
  - when a Frobenius helper really takes only `&curve`, prefer a dedicated
    extension trait method over a free function, rather than inflating broad
    core traits such as `EnumerableCurveModel`
  - when several finite-field group-order routes coexist for the same curve
    family, prefer one curve-side method such as `count_points(...)` with an
    explicit strategy enum, while keeping the underlying route visible in the
    returned report
  - when adding a Mestre-style route over `F_p`, prefer treating it as one
    more explicit `GroupOrderStrategy` with a report that preserves which side
    of Mestre's theorem won, rather than as a free-standing helper that hides
    the quadratic-twist branch
  - when that route needs small enum-conversion or dispatch helpers, keep
    those helpers crate-private unless they carry independent mathematical
    meaning for external callers
  - if a report already stores the step-by-step history of an iterative
    Frobenius-side algorithm, prefer deriving running summaries such as
    lower bounds from that history instead of storing duplicate cached totals
  - if one group-order strategy needs sampled points while the others are
    deterministic, prefer a sampler-aware sibling method such as
    `group_order_by_with_sampler(...)` over polluting the deterministic entry
    point with fake randomness defaults
  - if one internal step of a Frobenius-side algorithm can be improved from a
    naive search to BSGS without changing the surrounding mathematical report,
    prefer keeping the public report stable and treating the faster search as
    an internal implementation upgrade
  - if the Frobenius-side layer ends up with both naive and BSGS search
    engines for Hasse-interval multiples, prefer dedicated sibling modules for
    those engines and keep the trait surface as a thin caller-facing wrapper
  - if one internal BSGS implementation needs several orthogonal optimization
    knobs such as traversal order, fast negation, or parity information,
    prefer one crate-private config struct over a mutually exclusive enum
  - for those small internal config structs, prefer private fields and narrow
    `with_...` updaters over public struct literals at each call site
  - for experimental middle-out Hasse BSGS traversals, prefer separate
    left/right giant-step frontiers from the center block over alternating
    long re-centering jumps of one shared giant-step state
  - when evaluating a distribution-driven Hasse-search heuristic, prefer an
    adjacent deterministic corpus benchmark over many curve/point instances in
    the intended regime, rather than judging it only on one fixed instance
  - when comparing two internal search configurations for a Frobenius-side
    algorithm, prefer an ignored microbenchmark adjacent to that engine over a
    benchmark that pulls in unrelated curve layers
  - when a short-Weierstrass finite-field invariant can be phrased as
    `deg gcd(x^q - x, f(x)) > 0`, prefer computing `x^q mod f(x)` inside the
    cubic quotient ring over building the ambient degree-`q` polynomial
  - for finite-field cardinalities used in curve-side algorithms, prefer the
    validated ergonomic wrapper `FiniteField::order()` over repeating
    `cardinality().expect(...)`; and when building `H(q)`, prefer
    `HasseInterval::for_field::<F>()` or `HasseInterval::for_q(F::order())`
    so the call site stays honest about overflow and field validation
  - when a Hasse-interval BSGS search already knows the parity of `#E(F_q)`,
    prefer restricting the candidate progression to `M_0 + 2k` so the engine
    really pays for only one parity class, rather than checking the parity
    only after a full-width search
  - if that same invariant needs the cubic `x^3 + ax + b`, prefer exposing it
    as a method on `ShortWeierstrassCurve<F>` instead of rebuilding the
    coefficient vector ad hoc in each caller
  - if one Frobenius-side optimization prototype, such as a middle-out BSGS
    traversal, measures slower than the baseline, prefer keeping the baseline
    as default and documenting the prototype as explicit future work
  - when a group-order algorithm becomes mathematically central, such as the
    Mestre route, prefer adding property tests that compare it against the
    exhaustive `#E(F_p)` baseline on one prime field where both routes are
    still computationally honest
  - if a current short-Weierstrass algorithm implementation starts mixing
    field validation, twist selection, the main iteration, and final
    Frobenius/report packaging, prefer extracting small private helpers for
    those phases before introducing broader abstractions
  - if one Frobenius-side theory such as characteristic equations or isogeny
    invariance becomes a first-class namespace with its own reports and tests,
    prefer making it a direct `frobenius/<topic>/` submodule rather than
    wrapping it inside a second generic `verification/` layer
  - within `frobenius::torsion::matrix`, keep basis validation and subgroup
    coordinates separate from Frobenius-action/report assembly: coordinate
    tables belong with the torsion-basis side, while mod-`n` congruence
    checks belong with the final matrix-report construction
  - if that matrix submodule grows several small value types, prefer one file
    per mathematical object such as basis, matrix, and report, instead of a
    catch-all `types.rs`
  - if one of those congruence checks is really a statement made by a stored
    `FrobeniusTrace` package against another artifact, prefer making it a
    method on `FrobeniusTrace` rather than leaving it as a nearby free helper
  - when one route is mainly an implementation detail of that unified public
    method, prefer keeping the route-specific entry point crate-private rather
    than exposing two competing public doors
  - for barrel ergonomics, prefer keeping `elliptic_curves::...` and
    especially the crate root focused on stable curve-side values and primary
    entry points; route-specific step reports and diagnostics may stay public
    in their local module when needed, but they should not automatically be
    promoted into the broadest barrels
  - for the final Schoof stage, prefer keeping the CRT report as the source
    of truth for the prime-by-prime work and storing only the additional
    Hasse-resolution summary needed to decide whether the trace class is
    unique, ambiguous, or blocked
  - the natural automatic Schoof driver should own the stopping rule
    “accumulate odd primes until the CRT modulus exceeds the Hasse trace
    diameter bound `2⌊2√q⌋`”; keep the manually supplied odd-prime list only
    as an educational inspection surface
  - when a route-specific detailed report is generic in the field but the
    shared `GroupOrderReport` enum is intentionally non-generic, prefer one
    non-generic summary wrapper for the integrated strategy and keep the
    fully detailed generic report under its native Frobenius/Schoof namespace
  - if a route-specific report step type is mainly an internal algorithmic
    trace, prefer keeping that step type crate-private and exposing stable
    public summaries on the outer report such as labels, counts, or resolved
    candidates instead of raw step slices
  - arithmetic-derived educational views under `endomorphisms`, such as local
    conductor views, candidate local volcanic levels, or candidate-order cover
    relations, may be public when they are mathematically honest standalone
    surfaces; but downstream public reports in `isogenies` should still prefer
    stable summaries such as counts or possible levels over automatically
    re-exporting those fine-grained detail types
  - within `models::short_weierstrass`, prefer an austere root barrel:
    expose `ShortWeierstrassCurve` plus explicit submodules such as
    `point_order`, `group_order`, `group_exponent`, `division_polynomials`,
    `function_fields`, `isogenies`, and `isomorphisms`, rather than flattening
    their strategies and reports back into `short_weierstrass::*`
  - if `models::short_weierstrass::curve` grows mixed stories, prefer a
    directory-style module that separates the curve definition, invariants,
    twists/scalings, Frobenius actions, and tests, instead of one large
    `curve.rs`
  - if `models::short_weierstrass::curve/tests.rs` stops being tiny, prefer a
    `curve/tests/` directory split by mathematical story such as
    construction/invariants, isomorphisms, point enumeration, group law,
    order/structure, properties, and benchmarks, with shared fixtures in one
    local `shared.rs`
  - for shared short-Weierstrass secant/tangent formulas across backends,
    prefer a small internal `group_law_core/` module split into:
    point shape,
    backend ops trait,
    pure formula helpers,
    and one runner that owns `a` together with the chosen coordinate backend,
    rather than one flat file of free functions
  - within that same `group_law_core/` tree, prefer keeping slope and result
    reconstruction helpers in a `formulas/` submodule, while geometric control
    flow predicates such as vertical-opposite detection stay private to the
    runner that actually decides which branch of the group law to take
  - within `models::short_weierstrass::schoof::quotient_ring`, prefer the
    reduced quotient context, reduced quotient value, and partial-inverse /
    non-unit witness logic as separate sibling files rather than one growing
    mixed file; keep the barrel small and keep tests in a local `tests.rs`
  - within current `models::short_weierstrass::schoof` work, prefer keeping
    the reduced `(a(x), b(x) y)` endomorphism representation in its own file
    sibling to `quotient_ring/`, and keep polynomial substitution logic as a
    named helper on that value type rather than inlining Horner-style loops at
    each composition site
  - in the current Schoof additive arithmetic, remember that the additive zero
    endomorphism `P ↦ O` is not representable by one affine pair
    `(a(x), b(x) y)`; additive-combination helpers should therefore expose a
    dedicated zero branch instead of trying to encode it inside
    `ReducedEndomorphism`
  - in that same Schoof additive layer, detect additive-zero geometry before
    attempting quotient-ring inversion: both opposite affine maps
    `(a, by) + (a, -by)` and tangent-vertical doubling on the `y = 0` branch
    should return the zero endomorphism directly rather than surfacing a
    spurious non-unit denominator
  - for the current odd-prime Schoof driver before factor refinement exists,
    prefer one honest report that records `ψ_ℓ`, the reduced Frobenius data,
    the tested candidate residues, and whether the loop ended in
    `TraceFound`, `NonUnitDenominator`, or candidate exhaustion; do not hide
    the non-unit branch behind a generic error
  - for parity-sensitive division-polynomial recursions, especially the odd
    `ψ_{2m+1}` split by the parity of `m`, prefer checking one nontrivial
    output such as `ψ_5` or `ψ_7` against Sage or another exact CAS before
    building higher-level Schoof arithmetic on top of it
  - when that odd-prime Schoof report needs to expose reduced endomorphism
    data publicly, prefer reusing the canonical reduced arithmetic types
    `ReducedEndomorphism` and `ReducedEndomorphismAdditiveResult` rather than
    inventing report-local snapshot or mirror enums
  - for current Schoof quotient arithmetic, treat non-invertibility modulo the
    active polynomial `g(x)` as first-class algorithmic data: prefer an honest
    sum type carrying a gcd witness over collapsing that branch into a generic
    error
  - within `short_weierstrass::isogenies::velu::dual`, prefer moving tiny-field
    exhaustive witness searches such as “all scaling isomorphisms from this
    source curve to that target curve” onto `ShortWeierstrassCurve<F>` itself,
    and keep search-local duality checks as private methods on the Vélu search
    owner rather than as free helper functions
  - likewise, when `point_order` or `group_exponent` accumulate both API and
    several report/strategy types, prefer a thin `mod.rs` plus focused sibling
    files such as `api.rs`, `strategies.rs`, `reports.rs`, or
    `verification.rs`
  - within `point_order::from_multiple`, prefer splitting the story into
    report types, public API, validation/factorization bridge helpers, and the
    prime-by-prime recovery engine; if the final report can derive a quantity
    such as the remaining multiple from the exact order and recorded steps,
    prefer a derived accessor over storing a duplicate field
  - for `point_order` test-only support such as baselines or fixture helpers,
    prefer a sibling `tests_support.rs` or `tests/` tree over leaving that
    support nested under the production `from_multiple/` implementation module
  - within `models::short_weierstrass::group_exponent`, prefer separating the
    public strategy enum, random-point accumulation report, unified final
    report, and group-order-side verification into distinct sibling files, and
    keep verification-oriented curve methods in `verification.rs` rather than
    mixing them into the exponent-recovery API file
  - within `models::short_weierstrass::group_order`, prefer keeping
    deterministic curve-side dispatch in `api.rs`, the quadratic-character
    route in its own sibling file, and Mestre split by setup / loop / finalize
    responsibilities under `mestre/`, with sampler-aware docs explaining why
    `*_with_sampler(...)` exists instead of inventing hidden randomness
  - for analytic and Frobenius-heavy functionality in particular, prefer
    making callers spell the submodule (`analytic`, `frobenius`, etc.) once
    rather than flattening those specialized namespaces into the crate root
  - the same applies at the `elliptic_curves::` level itself: do not move
    order/group-exponent/group-order algorithms back into
    `models::short_weierstrass` just because the current public methods are
    implemented on `ShortWeierstrassCurve<F>`
  - exception: when an arithmetic subtree has in practice become entirely
    short-Weierstrass-specific, prefer moving its concrete reports and impl
    blocks under `models::short_weierstrass` and leaving only genuinely
    model-agnostic additive-group helpers under a small crate-private helper
    file near the `elliptic_curves` root, rather than re-introducing a broad
    `arithmetic/` helper subtree just for them
  - likewise, even core curve types such as `ShortWeierstrassCurve` should
    normally be reached through `elliptic_curves::...`, not through a
    duplicate crate-root alias
  - if one counting route is currently justified only for short-Weierstrass
    curves, prefer keeping that algorithm's implementation in the
    `short_weierstrass/` subtree until another concrete family really needs
    the same executable logic
  - when the backend already represents `F_{p^n}`, prefer reducing absolute
    Frobenius iterates modulo `n` in implementations and docs, since
    `π_p^n = id` on the represented field
  - when the backend already represents `F_q`, prefer a point-level helper
    named for the single relative map `π_q` rather than a powered helper with
    an observationally irrelevant iterate parameter
  - if multiple consumers need the same point-level relative Frobenius action,
    prefer a narrow capability trait such as `RelativeFrobeniusCurveModel`
    over hard-coding every helper to one concrete curve type
  - if exhaustive point counting is what enables a Frobenius helper, prefer a
    second stronger capability trait such as `FrobeniusTraceCurveModel`
    instead of folding enumeration requirements into the narrower point-level
    relative-Frobenius trait
  - when trace data is attached to a finite base field, prefer storing or
    exposing explicit `F_q` metadata such as `FiniteFieldDescriptor` over
    encoding only the characteristic `p`
  - when exposing the characteristic polynomial of the relative Frobenius,
    treat it as data derived from `FrobeniusTrace` with the convention
    `χ_{π_q}(T) = T^2 - tT + q`
  - when exposing the local zeta function, prefer deriving it from
    `FrobeniusCharacteristicPolynomial` rather than storing `q` and `t`
    again independently, so
    `Z(E/F_q, T) = (1 - tT + qT^2) / ((1 - T)(1 - qT))`
    keeps the characteristic polynomial as its single source of truth
  - when exposing the quadratic-twist point-count relation
    `#E(F_q) + #E'(F_q) = 2q + 2`, prefer a Frobenius-side report derived from
    the two `FrobeniusTrace` packages, while keeping the actual twist
    construction itself under `isomorphisms`
  - when exposing the fact that isogenous curves over `F_q` have the same
    point count and the same Frobenius trace, prefer a Frobenius-side report
    derived from the domain and codomain `FrobeniusTrace` packages of one
    explicit isogeny, rather than inflating the core `Isogeny` trait with
    extra arithmetic metadata
  - when extending that same experiment to an isogeny graph, prefer a report
    over the stored node representatives of the graph, compared against one
    reference node, rather than duplicating the graph layer's own edge-map
    reconstruction logic inside Frobenius helpers
  - when exposing counts over extensions `F_{q^n}` derived from Frobenius,
    prefer methods on `FrobeniusTrace` rather than on the curve directly, so
    the API keeps the “count once over `F_q`, derive many consequences later”
    story explicit
  - for that extension-count layer, prefer the `q`-relative recurrence
    `s_0 = 2`, `s_1 = t`, `s_n = t s_{n-1} - q s_{n-2}` with
    `#E(F_{q^n}) = q^n + 1 - s_n`, rather than reverting to prime-field-only
    notation in public docs or type names
  - for those derived extension counts, prefer arbitrary-precision integers
    over fixed-width `u128` / `i128`, so the educational API does not impose
    an artificial overflow boundary on `#E(F_{q^n})`
  - when a small represented extension field is still enumerable, keep the
    distinction explicit between:
    exhaustive point counting by `EnumerableCurveModel::order()`
    versus
    extension counts derived from a previously computed Frobenius trace
  - if both routes appear in one helper or test, prefer a named comparison
    report or clearly labeled variables so readers can see that one path is a
    direct enumeration witness and the other is the faster Frobenius-derived
    consequence being checked against it
  - when exposing Hasse-bound checks, prefer deriving them from
    `FrobeniusTrace` rather than duplicating `q`, `#E(F_q)`, and `t` as
    independent public state, and prefer the exact integral form
    `t^2 <= 4q` over floating-point approximations to `2 sqrt(q)`
  - when exposing finite-field point counts through
    `#E(F_q) = q + 1 + \sum_x χ(f(x))`, keep that route explicit as a
    character-sum algorithm distinct from the fully exhaustive point count
  - when exposing the search interval `H(q)` itself, keep it in the
    `frobenius/` layer next to `FrobeniusTrace` and compute its discrete
    integer endpoints with exact arithmetic, not floating-point square roots
  - when exposing a search for an annihilating multiple inside `H(q)`, prefer
    a Frobenius-side report and a curve-side method attached through the
    Frobenius capability layer, rather than treating that search as just one
    more point-count strategy
  - if exhaustive order search, order-from-multiple, and Hasse-driven order
    search all coexist for the same family, prefer one curve-side
    `point_order_by(...)` entry point with an explicit order strategy enum and
    route-preserving reports
  - if exhaustive and random-point routes coexist for the group exponent
    `λ(E(F_q))`, prefer one curve-side `group_exponent_by(...)` entry point
    with an explicit exponent strategy enum and route-preserving reports
  - if the Hasse-driven order route depends on how `#E(F_q)` was obtained,
    prefer carrying that choice explicitly in the order-strategy payload, for
    example through a `GroupOrderStrategy`, instead of hardwiring a hidden
    counting default
  - if the random-point exponent route depends on how sampled point orders are
    computed, carry that dependency explicitly in the exponent-strategy
    payload through a `PointOrderStrategy` instead of hardwiring one hidden
    order-recovery default
  - for route-preserving reports in this subtree, prefer deriving simple
    summary metadata such as counts or completion flags from one canonical
    stored payload when that keeps the report smaller and avoids duplicated
    state
  - if a point-sampling route only accumulates a lower bound for `λ(E(F_q))`,
    keep that accumulation separate from any later Hasse- or group-order-side
    verification of `#E(F_q))`; do not silently promote a lower-bound report
    into an exact-exponent report
  - likewise, if one route only certifies `#E(F_q)` by checking uniqueness in
    `H(q)` from another route's report, prefer keeping that certification as a
    separate verification helper rather than advertising it as a primary
    `GroupOrderStrategy`
  - for milestone examples that compare several group-order routes, prefer
    keeping them on one concrete curve over one prime field so the reader can
    compare exhaustive counts, character sums, Hasse searches, exponent
    witnesses, and Mestre side by side
  - when exposing ordinary versus supersingular classification, prefer
    deriving it from `FrobeniusTrace` via the general criterion `p | t`,
    where `p` is the base-field characteristic and `t` is the trace of `π_q`
  - in docs for that classification, it is good to mention the prime-field
    specialization `t = 0` for `F_p` with `p >= 5`, but keep `p | t` as the
    primary API-facing criterion so the same surface remains correct over
    extensions `F_{p^n}`
  - when first connecting Frobenius with torsion, it is acceptable to expose
    a report for the relative Frobenius `π_q` on exact-`n` rational torsion
    points, even though that report is tautological on `E(F_q)` because every
    listed point is already fixed by `π_q`
  - in docs for that current torsion/Frobenius layer, say explicitly that its
    value is pedagogical and preparatory for a later absolute-Frobenius view
    over extension fields, rather than pretending the first report is already
    a nontrivial arithmetic test
  - once that absolute-Frobenius view is added, keep the scope explicit:
    the first nontrivial helper may still be specialized to
    `ShortWeierstrassCurve<F>` if the codebase does not yet expose a generic
    absolute-Frobenius curve trait
  - if a Frobenius-on-torsion layer later exposes a matrix of `π_q` on a
    chosen basis of `E[n]`, keep that surface under `frobenius/`, not under
    `endomorphisms/`: it is a representation-on-torsion story, not an
    order/discriminant story
  - for that same matrix layer, keep the basis dependence explicit in both the
    type and the docs, and say directly that the matrix changes with the basis
    even though its trace and determinant modulo `n` do not
  - when the represented curve lives over a larger field than the trace base
    field, it is acceptable to realize `π_q` by one absolute `p^f`-power
    coordinate Frobenius, provided the docs also say explicitly that this
    requires the curve to be fixed by that Frobenius power
  - for base-defined curves viewed over `F_{p^r}`, it is good pedagogy to show
    both behaviors: torsion points fixed by `π_p` and torsion points moved by
    `π_p` but fixed by `π_p^r`
  - when a Frobenius/torsion report already records fixed-versus-moved points,
    it is good pedagogy to expose explicit derived accessors such as
    “already descends to `F_p`” versus “visible only in the chosen extension”
    instead of forcing readers to reinterpret raw fixed/moved counts on their own
  - when a Frobenius/torsion report already stores the point-to-image action,
    prefer exposing orbit decompositions directly from that stored data rather
    than making callers recompute Frobenius orbits from the curve a second time
  - when those stored point-to-image actions or orbit partitions are keyed by
    concrete finite-field points such as `AffinePoint<F>`, prefer small
    internal hashed lookup helpers over repeated linear scans, but keep that
    optimization behind internal APIs so broad public curve traits do not pick
    up unnecessary `Hash` bounds prematurely
  - once absolute Frobenius on torsion is available over `F_{p^r}`, it is also
    good pedagogy to expose the minimal positive fixing power `d` with
    `π_p^d(P) = P`, since that is the first direct bridge to “defined over
    `F_p`”, “defined over `F_{p^2}`”, and so on
  - Frobenius orbit helpers belong under `frobenius/` as their own small
    module, and should keep absolute versus relative orbit families explicit
    in both names and docs
  - for relative Frobenius orbits on `E(F_q)`, say directly that the current
    represented-field orbit story is trivial (singleton orbits)
  - for absolute Frobenius orbits on curves viewed over `F_{p^r}`, prefer
    exact orbit closure from the finite-field order bound over caller-supplied
    step caps when that bound is already known from the represented field
  - a standalone `FrobeniusDiscriminant` abstraction is now acceptable when
    it is derived directly from `FrobeniusTrace` and uses the shared
    integral `QuadraticDiscriminant` layer rather than duplicating ad hoc
    sign or squarefreeness logic
  - when exposing that Frobenius discriminant, prefer storing the originating
    `FrobeniusTrace` package plus the derived quadratic discriminant, rather
    than duplicating `q`, `#E(F_q)`, and `t` as independent public state
  - if a Frobenius-side helper returns the order `Z[π]`, name it explicitly as
    a Frobenius-generated order and say directly that it computes a natural
    suborder of `End(E)`, not the whole endomorphism ring
  - if a later report packages the sandwich
    `Z[π] ⊆ End(E) ⊆ O_K`, document it explicitly as a Frobenius-compatible
    candidate report and say directly that it does not identify the actual
    ring `End(E)` unless some later algorithm rules out the intermediate
    candidates
  - if that same report distinguishes ordinary from supersingular curves,
    keep the branches separate in both types and docs:
    the ordinary branch may use the imaginary-quadratic-order pipeline, while
    the supersingular branch should say directly that this module is not yet
    modeling the quaternionic endomorphism algebra as an order object
  - for branchy report enums in this endomorphism layer, prefer explicit
    variant names that encode the mathematical meaning, even if they are a bit
    longer, over short generic tags like just `Ordinary` or `Supersingular`
  - for small Frobenius-side reports, avoid duplicating synonymous accessors
    for the same derived payload; prefer one mathematically explicit name over
    parallel aliases like a generic `value()` plus a domain-specific accessor
  - when exposing an `\ell`-local endomorphism-ring view from
    `\Delta_\pi = v^2 D_K`, treat the global conductor `v` as the source of
    truth and expose the local datum as the valuation `v_\ell(v)`, not as a
    separate duplicated conductor object
  - in small endomorphism rustdocs where the notation is stable and common,
    prefer direct Unicode mathematical symbols such as `ℤ`, `π`, `Δ`, and `ℓ`
    over LaTeX-style command spellings inside code quotes
  - if we attach volcanic language to candidate orders, keep that layer
    explicitly arithmetic at first: `level_ℓ(O_f) = v_ℓ(f)` is a good derived
    value object under `endomorphisms`, but do not present it as a certified
    geometric level of the curve in an actual `ℓ`-isogeny volcano
  - that `\ell`-local view should live under `endomorphisms` as a derived
    Frobenius-side value object, while any generic integer-valuation helper
    used to compute it belongs under `numerics`
  - in docs for that local view, say explicitly that the global candidate
    orders form a divisibility poset, but after fixing one prime `\ell` the
    local exponents `0 <= b <= v_\ell(v)` do form a chain
  - when checking the characteristic equation pointwise, use the same
    relative-Frobenius convention
    `π_q^2(P) - [t]π_q(P) + [q]P = O`
    and do not silently replace `q` by `p` except in explanations of the
    prime-field special case
  - when a Frobenius helper is algorithmic rather than purely declarative,
    document its asymptotic cost explicitly in the rustdocs, especially when
    the current implementation is exhaustive over `E(F_q)` or iterates field
    Frobenius powers naively
  - for those asymptotic notes, prefer `Θ(...)` notation over `O(...)`, and
    say what the cost is measuring when that is not obvious from context, for
    example group additions/doublings, point-enumeration passes, or repeated
    field Frobenius updates
  - keep the distinction explicit between “lands on the Frobenius twist” and
    “is an endomorphism of the current curve”
  - use small extension-field examples to show that `E(F_p)` and `E(F_q)` can
    be distinguished by Frobenius fixed-point behavior even when the ambient
    coordinate backend is the same finite field
  - for Frobenius property tests, prefer algebraic invariants checked against
    small exhaustive witnesses, such as:
    trace-derived extension counts versus `curve.order()`,
    twist relations versus two direct point counts,
    isogeny invariance versus domain/codomain trace packages,
    or absolute-orbit partitions versus full enumerated point sets,
    instead of spending most of the property budget on `pretty()`/`Display`
    surfaces
- The `endomorphisms/` subtree is now an implemented part of the
  `elliptic_curves` architecture. For the current milestone, it is acceptable
  and preferable to:
  - keep the scope finite-field-first and short-Weierstrass-first if that
    matches the already implemented Frobenius infrastructure

- The short-Weierstrass `function_fields/` subtree under
  `models::short_weierstrass` is now also an implemented part of the
  `elliptic_curves` architecture. For the current milestone, it is acceptable
  and preferable to:
  - keep the scope explicitly short-Weierstrass-first
  - model `F(E)` via the basis `1, y` over `F(x)` as pairs `(A(x), B(x))`
  - document directly that multiplication uses the specific relation
    `y^2 = x^3 + ax + b`
  - keep short-Weierstrass-specific substitution helpers such as evaluating
    `x^3 + ax + b` at a function-field element under
    `elliptic_curves::short_weierstrass::function_fields`, not in downstream
    isogeny helpers
  - keep the ambient curve as explicit runtime context on both the family and
    the stored function values instead of pretending the current `Field` trait
    can express a runtime-dependent curve family
  - when the current algebra really needs runtime ambient data, prefer
    implementing a dedicated ambient trait such as
    `fields::traits::AmbientField`
    rather than pretending the family is type-level static
  - use the conjugation and norm formulas pedagogically, including the inverse
    formula `(A, B)^(-1) = (A / (A^2 - fB^2), -B / (A^2 - fB^2))` whenever the
    norm is non-zero
  - build educational value objects or reports from existing data such as
    `FrobeniusTrace`, `FrobeniusCharacteristicPolynomial`, and the Frobenius
    discriminant `t^2 - 4q`
  - say explicitly when a result is only reporting the quadratic order
    suggested by the current Frobenius data, rather than claiming a proved
    computation of the full geometric endomorphism ring or its maximal order
  - keep the distinction explicit between:
    Frobenius as one concrete endomorphism over finite fields
    versus
    the whole ring `End(E)` or `End(E) ⊗ ℚ`
  - prefer narrow names such as “order”, “candidate order”, “discriminant”,
    or “Frobenius-generated suborder” when those are the honest mathematical
    objects currently supported, instead of using stronger names prematurely
  - keep graph heuristics, volcano language, and endomorphism-ring claims
    separate unless the implementation really computes the arithmetic object
    being named
  - document algorithmic cost in rustdocs with `Θ(...)` notation just as in
    the existing Frobenius layer, and say clearly whether the cost is only
    integer arithmetic on stored trace data or still depends on exhaustive
    curve enumeration upstream
  - for the first arithmetic helper under `endomorphisms/`, prefer a small
    integer-backed `QuadraticDiscriminant` value object before any fuller
    quadratic-order or quadratic-field abstraction
  - that first discriminant layer should focus on lightweight integral
    arithmetic and classification facts such as sign, congruence modulo `4`,
    squarefree/fundamental status, and construction from Frobenius data
    `D = t^2 - 4q`, without pretending it already models a full order or the
    whole CM field
  - if that discriminant layer needs a shared squarefreeness predicate or
    other domain-agnostic integer arithmetic, prefer moving the helper to
    `numerics/` rather than duplicating it locally under `endomorphisms/`
  - if `endomorphisms/` grows past one source file, prefer a dedicated
    `tests.rs` for the subtree's shared unit tests instead of scattering the
    first arithmetic checks across multiple tiny inline test blocks
  - if that same discriminant layer grows to the canonical decomposition
    `Δ = v^2 D_K`, prefer a dedicated value object for the factorization that
    stores the original discriminant, the positive square root factor `v`, and
    the resulting fundamental discriminant `D_K` explicitly, instead of
    returning a loose tuple
  - once that factorization is used to enumerate all intermediate quadratic
    orders between `ℤ[π]` and `O_K`, prefer a dedicated candidate-set value
    object whose public story is “divisors of the conductor `v`”, with the
    orders returned in one deterministic order rather than as an unordered bag;
    document explicitly that this list is only a storage/view convention, while
    the mathematically natural structure is the divisibility poset on the
    conductors, equivalently the order-containment poset on the corresponding
    quadratic orders
  - if the canonical factorization logic grows branchy, prefer private
    helper extraction by mathematical case such as the `Δ ≡ 1 mod 4` and
    `Δ ≡ 0 mod 4` branches, so the public factorization entrypoint reads as a
    short dispatcher rather than one long mixed-case implementation
  - if the next layer models an imaginary quadratic order
    `O_f = Z + f O_K`, prefer constructing it explicitly from the pair
    `(D_K, f)` or from the existing canonical factorization `Δ = f^2 D_K`,
    and keep the pipeline
    `QuadraticDiscriminant -> QuadraticDiscriminantFactorization -> ImaginaryQuadraticOrder`
    visible in the public API rather than hiding it behind one oversized
    constructor
  - for the first radicand-based quadratic-field bridge, keep the scope
    intentionally narrow: start only from integer radicands `m < 0`, reduce
    to the squarefree part, apply the classical mod-`4` rule, and expose only
    the resulting fundamental discriminant `D_K` and maximal order `O_K`
    publicly; do not introduce a broader public quadratic-field or
    integral-basis API until a concrete consumer needs it
  - when exposing relative indices between two such orders in the same field,
    prefer the direct conductor formula
    `[O_{f_2} : O_{f_1}] = f_1 / f_2` under the inclusion
    `O_{f_1} ⊆ O_{f_2}`, and keep any first API surface small, for example a
    method on `ImaginaryQuadraticOrder`, before introducing a larger report type
  - if a candidate-set helper is added on top of that, prefer a very small
    convenience method specialized to indices over `Z[π]`, reusing the
    general order-index method rather than duplicating the arithmetic
  - if a pedagogical poset or Hasse-diagram view is added for candidate
    orders, prefer storing the mathematical cover relations and edge labels in
    `endomorphisms/`, while keeping the human-readable textual rendering under
    `visualization/elliptic_curves`
    constructor
  - when that same factorization is also used to recover the maximal order
    `O_K`, prefer exposing a dedicated helper such as `maximal_order()` on the
    factorization object itself, so the relationship between `Z[π]` and
    `O_K` stays visible through the shared `D_K`
- Complex-analytic scaffolding may introduce small numerical
  helper types when the docs stay explicit that the current goal is
  educational floating-point experimentation rather than numerically certified
  complex analysis.
  Small modular value objects such as a stored `q = e^{2π i τ}` parameter are
  acceptable when they keep the upper-half-plane input explicit and make later
  `q`-expansion code easier to explain.
  For first `q`-expansion approximations, it is acceptable to ship only a
  short explicit coefficient table when the docs say clearly which terms are
  currently implemented and how truncation counts those terms.

## Design priorities

- Mathematical honesty before feature count.
- Clear point representations before group-law optimization.
- Conservative public APIs that explain their preconditions.
- Small, verifiable steps.

## Representation rules

- Prefer representations that make invalid states hard to express.
- The point at infinity should be modeled explicitly rather than smuggled
  through meaningless affine coordinates.
- Representation-level coordinate transports that do not themselves certify
  curve membership belong on `AffinePoint` rather than being hidden inside one
  feature module. Document clearly that they move coordinates only and that
  target-curve validation remains the caller's responsibility.
- If a constructor claims to return a point on a curve, it should validate that
  claim.
- If a curve model is only valid away from special characteristics, document
  that fact directly in the type or constructor docs.

## Scope guidance

- It is fine to start with affine membership checks, discriminants, and simple
  point constructors.
- Model-specific invariants can stay as inherent methods when they belong only
  to one presentation, such as short-Weierstrass invariants on
  `ShortWeierstrassCurve`.
- When an invariant is useful across multiple consumers but still not universal,
  prefer a narrow capability trait over adding it directly to `CurveModel`.
- Short-Weierstrass coefficient-scaling helpers such as `scaled_by` and
  `isomorphic_via_scale` belong on `ShortWeierstrassCurve` itself, since they
  describe model-specific coefficient transport rather than a generic
  curve-morphism interface.
- Explicit short-Weierstrass scaling witnesses belong in
  `ShortWeierstrassIsomorphism<F>` and currently cache their codomain curve as
  validated state. That means the `CurveIsomorphism` trait may return the
  codomain by reference honestly, without recomputing it.
- For short-Weierstrass isomorphisms, use the convention
  `\phi_u : E -> E'`, `(x, y) -> (u^2 x, u^3 y)`.
  If `E : y^2 = x^3 + ax + b`, then document and implement the image model as
  `E' : y^2 = x^3 + a'x + b'` with `a' = u^4 a` and `b' = u^6 b`.
  Treat this as the canonical normalization for the current isomorphism layer unless a later extension
  explicitly introduces a second convention and explains the translation.
- Keep the distinction explicit between:
  same `j`-invariant = isomorphic over an algebraic closure,
  versus
  isomorphic over the base field = there exists `u in F^*` with
  `a' = u^4 a` and `b' = u^6 b`.
  Do not collapse those two notions in docs or API names.
- For quadratic-twist objects, store the twist factor `d` as the primary data.
  Do not store a base-field scaling witness `u` as mandatory state, since that
  witness may fail to exist exactly in the genuinely quadratic case.
- For Frobenius helpers, keep the distinction explicit between:
  coefficients fixed by the prime-field Frobenius,
  versus
  points fixed only after the full relative Frobenius of the chosen finite
  base field.
  Do not collapse “the curve is defined over `F_p`” with “the point is
  `F_p`-rational” in docs or API names.
- Point enumeration is acceptable only when the base field is explicitly small
  and enumerable. Say so in docs.
  This includes short-Weierstrass curves over small enumerable extension
  fields once the field backend honestly provides both full element
  enumeration and square-root discovery.
- Division-polynomial torsion search is acceptable only when the field
  capabilities are honest about what is being used:
  - `EnumerableFiniteField` for exhaustive `x` or point scans
  - `SqrtField` only when the code genuinely uses square-root lifting
- For `EnumerableFiniteField`, exhaustive witness search is acceptable for
  pedagogical helpers such as finding a concrete short-Weierstrass scaling
  isomorphism or enumerating all compatible base-field isomorphisms between two
  short-Weierstrass models. Keep the docs honest that this is a small-field
  educational routine, not a large-field optimized algorithm.
- Point-order and torsion helpers are acceptable only when the ambient group is
  explicitly small and enumerable. Say directly that the current algorithms
  use direct traversal or repeated addition rather than efficient large-group
  techniques.
- Generic exact-order helpers such as `point_has_exact_order(...)` and
  `points_of_exact_order(...)` should live as default methods on
  `FiniteGroupCurveModel`, not in a generic shared bucket or a free-standing
  torsion module.
- Division-polynomial-based torsion helpers belong in
  `src/elliptic_curves/models/short_weierstrass/division_polynomials/`, and
  should keep the distinction explicit between:
  - `x`-candidates
  - torsion candidates from `ψ_n(P)=0`
  - torsion points after extra validation
  - exact-order-`n` torsion points
- Keep odd/even division-polynomial semantics explicit in both naming and
  docs:
  - odd `ψ_n` live directly in `F[x]`
  - even `ψ_n` have the shape `y ε_n(x)`
  - for even `n`, `ψ_n(P)=0` means `y(P)=0` or `ε_n(x(P))=0`
- When callers need only an `x`-coordinate test, prefer one explicit helper
  that dispatches by parity to either `ψ_n(x)` or the stripped even factor
  `ε_n(x)` instead of open-coding the odd/even branch in multiple places.
- For even division polynomials, keep the documentation explicit that
  `ψ_n = y ε_n(x)` and that the `y = 0` branch can contribute lower-order
  torsion candidates.
- Do not rush into optimized formulas, scalar multiplication, serialization, or
  cryptographic hardening.
- If a model exposes a group-law trait, keep the docs explicit about which
  operations are baseline educational formulas and which errors represent
  invalid off-curve inputs.
- If a new curve API depends on extra field capability, such as square roots,
  prefer a narrow trait bound like `SqrtField` over broadening unrelated base
  traits.
- For analytic floating-point helpers, prefer small explicit value objects
  such as tolerances or normalization settings over hidden global constants.
  Document what the knobs mean and keep preset constructors easy to compare.
  When those helpers are shared with `fields` or other domains, prefer placing
  them in sibling numerical infrastructure and reexporting them here instead of
  making `fields` depend on `elliptic_curves`.
  For upper-half-plane and similar analytic domain types, document explicitly
  how near-boundary floating-point inputs are treated under the default
  tolerance policy.
  Named educational sample points such as `i`, `rho`, or a generic interior
  example are acceptable when they help later modular or lattice examples stay
  concrete and readable.
  For complex lattices, it is acceptable to require a positively oriented
  basis at construction time when that keeps the associated `tau` parameter
  honest and avoids hidden reordering conventions.
  For boxed lattice enumeration helpers, store both the integer indices and
  the concrete complex value so examples can show the arithmetic transparently.
  For lattice-sum truncation helpers, prefer a validated value object with a
  private stored radius rather than passing raw `usize` knobs everywhere.
  If the truncation models a square index box in `ℤ²`, say so directly in the
  docs and keep zero-radius rejection explicit when the origin-only box would
  be mathematically unhelpful for the intended analytic routine.
  For torus and fundamental-parallelogram helpers, it is acceptable to model
  a torus point canonically by reduced coordinates in a half-open region, as
  long as the docs explain clearly that the meaning is still relative to a
  chosen lattice.
  When a helper type claims to represent canonical coordinates in that region,
  prefer validating the constructor and keeping the stored coordinates behind
  explicit accessors instead of exposing unchecked public fields.
  When reducing floating-point coordinates modulo the unit square, document
  any boundary-snapping convention explicitly, especially if values very close
  to `0` or `1` are normalized back to `0`.
  When an analytic module grows a mixed set of responsibilities, prefer
  promoting it from a single `foo.rs` file to a `foo/` module with focused
  subfiles such as basis, coordinates, reduction, or reports, and keep a
  dedicated `tests.rs` when the test surface is no longer tiny.
  If equality is mathematically meaningful only relative to ambient lattice or
  curve context, prefer an explicit `..._eq(...)` method on that ambient type
  over deriving `PartialEq` on the context-free value alone.

## Error conventions

- Keep recoverable curve-domain failures in `CurveError`.
- Prefer specific variants such as unsupported characteristic, singular curve,
  or point-not-on-curve over ad hoc strings.
- Add a new error variant only when it expresses a genuinely distinct curve
  failure mode.
- For analytic lattice helpers, keep “degenerate basis” separate
  from “non-positive orientation” when the code or docs need to explain why a
  basis was rejected pedagogically.
- Torsion-order validation errors that are generic to curve-group logic, such
  as “order must be positive”, belong in `CurveError` rather than in a
  feature-local error enum.
- For analytic torus torsion, it is acceptable to expose the analytic
  counterpart of `E[n]` through reduced lattice indices `(a, b; n)` with
  explicit docs for `E[n] ≅ (1/n)Λ / Λ`, provided the constructor validates
  `n > 0` and `0 ≤ a, b < n`.
- When the torus-side torsion API distinguishes all `n`-torsion from
  primitive `n`-torsion, document directly that “primitive” means exact torus
  order `n`, equivalently `gcd(a, b, n) = 1`.
- If a torus torsion point stores the arithmetic index, reduced fundamental
  coordinate, and explicit complex representative together for pedagogy,
  prefer private fields plus accessors so those redundant views cannot drift
  out of sync.

## Testing expectations

- Test both valid and invalid curve construction.
- Test both valid and invalid point construction.
- Test the point at infinity behavior explicitly when it participates in the
  public model.
- When group operations are exposed, test identity, inverses, doubling, scalar
  multiplication, and at least one small exact associativity example.
- When point-order or torsion helpers are exposed, test at least one identity
  case, one non-trivial finite-order example, and one invalid off-curve input.
- For rational torsion over `Q`, prefer Unicode mathematical notation in
  rustdocs when it improves readability, such as `E(Q)_tors`, `ℤ/nℤ`, `ψₙ`,
  `Δ`, and `ℓ`. Each implementation pass in that area should update the
  relevant `AGENTS.md` guidance when it introduces a new local convention, and
  verification should stay module-focused first, for example
  `cargo test -q rational_torsion`, before broadening to neighboring
  short-Weierstrass tests only when the touched surface crosses those
  boundaries.
- For rational-torsion value objects that represent a certified mathematical
  classification, prefer one validating constructor over public enum variants
  that can encode impossible states. Keep any looser candidate shape as a
  separate input type, and let the validated report store only the checked
  value object. Once such a value object has been constructed, avoid exposing
  extra `has_valid_shape`-style predicates that re-check the invariant; reserve
  validation for the constructor boundary. Classification from exact verified
  point orders should live on `RationalTorsionGroup`, not in the verification
  pipeline that merely computes those orders.
- For rational-torsion reports, keep candidate accounting as a checked
  invariant at construction time. Store only the canonical point payload plus
  the total candidate count, and derive rejected counts from those values
  instead of storing a second mutable-looking aggregate.
- For rational-torsion verification over `Q`, avoid multiplying every
  non-torsion candidate by the full Mazur exponent `27720`: exact rational
  coordinates can grow too quickly. Prefer testing the finite list of possible
  non-identity Mazur point orders `2, ..., 10, 12`, with the identity handled
  separately as order `1`; this is equivalent for rational torsion
  classification and keeps focused tests small. Keep Mazur-specific constants
  in a small Mazur module rather than coupling classification or verification
  to Lutz-Nagell enumeration. Once a candidate's exact order has been computed
  during verification, carry that order forward into group classification
  instead of recomputing it. Do not add a separate defensive closure under
  negation to this route: the Lutz-Nagell candidate report already enumerates
  `y = 0` and both signs `±y` whenever `y² | Δ`, and it already deduplicates
  candidates. Prefer reusing the shared affine-point sorter for deterministic
  report ordering instead of exposing lower-level point comparators.
- For staged rational-torsion scaffolding, keep modules, constants, reports,
  candidate-shape enums, and helper witnesses `pub(crate)` or narrower until a
  tested curve-side entry point needs to expose them. Public API should grow
  from concrete `ShortWeierstrassCurve<Q>` workflows, not from placeholder
  internals. Once rational torsion uses a dedicated `tests/` module directory,
  keep tests split by responsibility, for example classification, reports,
  integral models, and errors, instead of letting `tests/mod.rs` become a
  catch-all.
- For rational-torsion examples, prefer the curve-side
  `ShortWeierstrassCurve<Q>::rational_torsion()` entry point over exposing
  integral-model or Lutz-Nagell internals. Keep the public report inspectable,
  but leave construction and verification helpers crate-private.
- For rational-torsion regression tests, include at least one case where
  Lutz-Nagell produces non-torsion candidates that are rejected, and include
  cyclic Mazur shapes beyond the early `2`-torsion and `ℤ/6ℤ` examples so the
  verification route is visibly not specialized to those cases.
- For rational-torsion integral-model witnesses, reuse the existing
  `ShortWeierstrassIsomorphism<Q>` / `scaled_by(u)` convention for
  `ϕᵤ : E → E'` instead of storing a second independent `(curve, scale)` pair.
  The scaled integral companion should be the isomorphism codomain, so the
  coordinate-change data remains the single source of truth. When computing
  the integral scale for `E: y² = x³ + ax + b`, choose `u` from denominator
  prime-power valuations via the shared `numerics` rational-denominator
  clearance helper, with exponent `4` for `a` and `6` for `b`; transport points
  through the stored isomorphism and its inverse rather than repeating
  coordinate formulas in the torsion layer.
- For the rational-torsion Lutz-Nagell route, keep the candidate report on the
  integral companion model and validate every returned affine candidate with
  `curve.point(...)`. Use the polynomial-side integer-root helper for
  equations such as `x³ + Ax + B - y² = 0` instead of embedding a local
  rational-root/divisor search in the curve module. When a deterministic point
  order is needed for reports or tests, prefer a direct comparator over
  allocating cloned coordinate sort keys. Document the candidate-enumeration
  cost in `Θ(...)` notation, naming reused subroutine costs such as integer
  factorization or polynomial-side integer-root search, but keep the expression
  coarse enough to read comfortably in rustdocs. Prefer constructing the
  candidate report through an associated constructor such as
  `LutzNagellCandidateReport::from_integral_model(...)` instead of exposing a
  parallel free function that returns the same report.
- For division-polynomial torsion helpers, test at least:
  - one explicit low-degree base formula
  - one recursive odd case and one recursive even case
  - one property check of the form `P` non-trivial `n`-torsion implies
    `ψ_n(P) = 0`
  - one comparison against exhaustive enumeration
  - one case that distinguishes raw candidates from exact-order points
  - one explicit negative case where an odd division polynomial does not
    vanish on a generic non-torsion point
  - one explicit root-lifting case and one explicit non-liftable-root case
- When a helper depends on field-side capabilities, add at least one test that
  exercises the positive path and one that shows the honest negative path.
- For analytic numerical helper types, test the preset constructors
  directly and keep the expected constants explicit in the tests.
- For enumeration helpers, test the identity case, finite-point count, and at
  least one small exact order example.
- For short-Weierstrass isomorphism comparisons, include explicit tests for
  the special `j = 0` (`a = 0`) and `j = 1728` (`b = 0`) families, in
  addition to generic `a,b != 0` examples.
- When short-Weierstrass isomorphisms are cached objects rather than
  recomputed witnesses, test both the coefficient transport and the cached
  codomain-facing behavior through the generic `CurveIsomorphism` trait
  surface.
- For short-Weierstrass automorphism helpers over enumerable finite fields,
  test the generic case separately from the `j = 1728` and `j = 0` special
  families, since those special loci can admit extra automorphisms.
- When a directory-style module already separates behavior by responsibility,
  it is acceptable to extract the shared structs into a local `types.rs` once
  `mod.rs` starts acting mainly as a type bucket plus submodule list.
- For quadratic-twist helpers, test at least:
  preservation of the `j`-invariant,
  one square-factor case that stays base-field isomorphic,
  one non-square case that is not base-field isomorphic in the chosen sample,
  and the point-count relation `#E(F_p) + #E^(d)(F_p) = 2p + 2` over small odd
  prime fields.

## Documentation expectations

- Public curve items should explain the mathematical model they represent.
- If a formula is valid only in characteristic different from `2` and `3`, say
  so directly.
- If a feature is educational, partial, or not yet a full group-law layer, say
  so explicitly.
- If an invariant is attached to a specific curve presentation, document both
  its defining formula and its mathematical role.
- If a helper only makes sense for small finite fields, say so directly in the
  rustdocs.
- Use concrete examples where they clarify the model.
- When documenting Frobenius helpers, explain both:
  - whether the map is absolute `π_p` or relative `π_q`
  - whether the output lies on the same curve or on a Frobenius twist
  - and prefer the `π_p` / `π_q` notation consistently once that notation is
    introduced in the surrounding text
- When documenting division-polynomial helpers, explain both:
  - the algebraic recurrence being implemented
  - the computational cost under the current dense naive multiplication backend
- When documenting analytic numerical helpers, state clearly whether a
  constructor validates its inputs or merely packages caller-supplied
  tolerances for later algorithms.

## Review heuristics

A good change under `src/elliptic_curves` should improve at least one of:

- invariant safety
- readability
- mathematical honesty
- test coverage

If a curve change makes the point model or equation semantics harder to
explain, it is probably moving too fast for the current phase.

- Within `elliptic_curves::frobenius`, avoid a leftover `maps/` namespace once
  concrete curve-side actions have moved onto `ShortWeierstrassCurve<F>`.
  Prefer placing pure metadata objects such as `AbsoluteFrobenius` and
  `RelativeFrobenius` directly under `frobenius/metadata.rs`, and orbit value
  types plus small crate-private orbit utilities directly under
  `frobenius/orbit.rs`.
- Frobenius-facing reports should expose exact integer data as
  `BigUint`/`BigInt` rather than fixed-width `u128`/`i128`. If a staged
  algorithm still needs a primitive value, keep that conversion behind a
  deliberately named compatibility helper and do not let the primitive type
  leak back into the public report surface.
- Frobenius and torsion tests over extension fields should use the smallest
  field that still exercises the intended phenomenon, and should cache
  enumerated point sets when several assertions share the same curve. Keep
  larger extension-field examples for ignored or integration-style coverage
  unless their extra size is mathematically necessary.
- Hasse-search engines may use local memory-sized indices for baby-step tables
  or benchmark counters, but their public interval, report, and returned
  annihilating-multiple surfaces should remain `BigUint`/`BigInt`.
- Within `elliptic_curves::frobenius::torsion`, prefer splitting the pointwise
  torsion-action story from the matrix-on-`E[n]` story, and when a public
  entry point is specifically a short-Weierstrass computation, prefer exposing
  it as a `ShortWeierstrassCurve<F>` method rather than as a free function in
  the torsion namespace.
- Within `elliptic_curves::frobenius::cm`, keep arithmetic candidate generation
  separate from curve-side certification. When a helper turns a CM candidate
  into curve-specific trace information, prefer a curve-side trait method over
  a free function and document the witness being used, such as a point `P`
  distinguishing `[p + 1 - t]P` from `[p + 1 + t]P`. For retrying witnesses,
  prefer explicit point iterators, small-field `EnumerableCurveModel`
  enumeration, or the existing `PointIndexSampler` abstraction rather than
  adding a direct randomness dependency.
- Within `models::short_weierstrass::curve`, if Frobenius helpers start mixing
  relative action, absolute twist transport, and orbit enumeration, prefer a
  `curve/frobenius_actions/` directory split by those stories rather than one
  flat file.
- Within `models::short_weierstrass::division_polynomials`, prefer asking the
  curve a higher-level question such as whether `ψ_n` carries a `y`-factor
  instead of exposing a public criterion-dispatch enum just to let downstream
  code branch on parity-shaped behavior.
- Within `models::short_weierstrass::function_fields`, prefer keeping the
  runtime ambient object under `field/`, the raw `A(x) + yB(x)` value story
  under `value/`, and the generic-point / function-field group-law story
  under `point/`, with `tests/` split along those same narratives rather than
  one monolithic file.
- In that same `function_fields` tree, if affine-point validation is needed in
  more than one constructor/helper, prefer one crate-private validation helper
  reused by both `affine(...)` and later consistency checks instead of
  duplicating the same curve/equation logic.
- For cyclic-group algorithms from educational problem sets, such as prime
  `r`-th root extraction in a finite cyclic group, keep the generic integer
  and route-report objects under `elliptic_curves::group_algorithms` first.
  Adapt concrete curve families to that shared core later instead of embedding
  the algorithm directly in one curve model.
- Within `group_algorithms::cyclic_roots`, keep setup decomposition, Bezout
  data, discrete-log search steps, route traces, outcomes, reports, and tests
  in separate files. The module root should stay as a small index rather than
  becoming the implementation body.
- In that same cyclic-root setup layer, derive dependent invariants such as
  `a`, `r^k`, and `k` from the canonical input pair `(|G|, r)` instead of
  accepting them independently from callers.
- Keep the executable prime-root route in `group_algorithms::cyclic_roots::algorithm`
  close to the exercise statement: compute `α = aγ`, `β = r^kγ`, brute-force
  `α = xδ`, and then apply `ρ = s(x/r)δ + tβ`. Use small cyclic curve groups
  as regression tests before introducing the large problem-set curve example.
- For that same route, include at least one regression test that checks the
  recorded trace against the exercise identities, not only the final predicate
  `[r]ρ = γ`; this keeps `α`, `β`, the brute-force `x`, and Bezout data from
  becoming decorative report fields.
- Expose that cyclic-root route to curve consumers through the
  `CyclicGroupPrimeRootCurveModel` capability trait and
  `cyclic_group_prime_root(...)` method. Treat the cyclic-group hypothesis and
  the claim that `δ` generates the full `r`-Sylow subgroup as external data or
  a separate certificate; the algorithm itself does not prove group cyclicity.
  The free function in `algorithm.rs` is an engine, not the primary API shape.
- Before adding the large problem-set curve example, keep the public
  cyclic-root surface limited to the curve method plus read-only report/value
  accessors, and protect that boundary with an integration test that imports
  `elliptic_curves::group_algorithms::cyclic_roots` from outside the crate.
- When adding the large Problem Set 2 cyclic-root example, keep all external
  constants visible in the example file, document which scalar `c` is being
  used, derive each `δ = (n/r)P` in code, and assert `[r]R = γ` for every
  printed root so the output is a checked computation rather than a transcript.
- Do not introduce a broad standalone group trait just to stage the first
  cyclic-root implementation. Start with the existing `GroupCurveModel`
  operations and extract a truly generic additive/multiplicative group action
  trait only when a second non-curve consumer such as `F_q^×` is ready.
