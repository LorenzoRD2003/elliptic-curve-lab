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
- The first function-field layer should stay short-Weierstrass-specific until
  the repo has a genuinely justified abstraction for other curve models.
- Public substitution helpers that are intrinsic to the short-Weierstrass
  function field, such as evaluating polynomials or rational functions in the
  distinguished `x`-coordinate at a function-field element, belong under
  `elliptic_curves::function_fields` rather than under downstream consumers
  such as `isogenies`.
- Likewise, generic-point arithmetic in `E(F(E))` belongs under
  `elliptic_curves::function_fields`: if Vélu, scalar multiplication, or other
  consumers need rational addition, translation, or doubling formulas on the
  generic point, prefer exposing one shared point representation and helpers
  there instead of re-encoding secant/tangent formulas in each consumer.
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
  - expose staged public APIs such as:
    `rational_x_candidates_for_division_polynomial(...)`,
    `torsion_candidates_from_division_polynomial(...)`,
    `torsion_points_from_division_polynomial(...)`,
    `exact_n_torsion_points_from_division_polynomial(...)`, and
    `compare_division_polynomial_torsion_with_enumeration(...)`
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
  - when several finite-field point-count routes coexist for the same curve
    family, prefer one curve-side method such as `count_points(...)` with an
    explicit strategy enum, while keeping the underlying route visible in the
    returned report
  - when one route is mainly an implementation detail of that unified public
    method, prefer keeping the route-specific entry point crate-private rather
    than exposing two competing public doors
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
- A future `endomorphisms/` subtree is now an intended part of the
  `elliptic_curves` architecture. For the first milestone, it is acceptable
  and preferable to:
  - keep the scope finite-field-first and short-Weierstrass-first if that
    matches the already implemented Frobenius infrastructure

- A `function_fields/` subtree is now also an intended part of the
  `elliptic_curves` architecture. For the current milestone, it is acceptable
  and preferable to:
  - keep the scope explicitly short-Weierstrass-first
  - model `F(E)` via the basis `1, y` over `F(x)` as pairs `(A(x), B(x))`
  - document directly that multiplication uses the specific relation
    `y^2 = x^3 + ax + b`
  - keep short-Weierstrass-specific substitution helpers such as evaluating
    `x^3 + ax + b` at a function-field element under
    `elliptic_curves::function_fields`, not in downstream isogeny helpers
  - keep the ambient curve as explicit runtime context on both the family and
    the stored function values instead of pretending the current `Field` trait
    can express a runtime-dependent curve family
  - when the current algebra really needs runtime ambient data, prefer
    implementing a dedicated ambient trait such as `fields::AmbientField`
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
  `points_of_exact_order(...)` belong in `src/elliptic_curves/torsion.rs`.
- Division-polynomial-based torsion helpers belong in
  `src/elliptic_curves/division_polynomials/`, and should keep the distinction
  explicit between:
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
