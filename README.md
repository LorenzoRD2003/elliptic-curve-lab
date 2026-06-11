# elliptic-algorithms-lab

An educational Rust library for finite fields, polynomials, elliptic curves,
and related algebraic or cryptographic algorithms.

The project is being built in public, in small steps, with an emphasis on:

- mathematical honesty
- readable Rust
- tests for algebraic behavior
- educational formatting and visualization helpers
- architectures that can grow without becoming confusing

Current state:

- the finite-field division-polynomial layer is complete enough to study educational division polynomials,
  rational torsion over small finite fields, and their interaction with
  Vélu-style finite-field workflows
- the first Frobenius layer now distinguishes:
  - absolute Frobenius `π_p` on coordinates
  - relative Frobenius `π_q` over a chosen finite base field
  - Frobenius twists of short-Weierstrass curves
  - a first Frobenius-trace layer based on exhaustive point counting and the
    identity `t = q + 1 - #E(F_q)` for small finite base fields
  - a Frobenius-discriminant report for
    `Δ_π = t^2 - 4q`, built directly from the same trace package
  - the local zeta function `Z(E/F_q, T)` derived explicitly from the same
    Frobenius characteristic polynomial
  - the quadratic-twist Frobenius relation
    `#E(F_q) + #E'(F_q) = 2q + 2` together with the equivalent trace
    negation `t' = -t`
  - the isogeny Frobenius relation showing that explicit isogenies preserve
    both `#E(F_q)` and the Frobenius trace
  - the corresponding graph-level Frobenius experiment showing that all stored
    node representatives in a small isogeny graph share the same trace
  - derived extension counts `#E(F_{q^n})` from the same trace via the
    classical order-2 Frobenius recurrence
  - exact Hasse-bound verification via the equivalent integer inequality
    `t^2 <= 4q`
  - ordinary versus supersingular classification derived from the Frobenius
    trace through the criterion `p | t`
  - a first Frobenius/torsion bridge for exact rational `n`-torsion under the
    relative Frobenius `π_q`
  - a first nontrivial Frobenius/torsion bridge for absolute Frobenius `π_p`
    on torsion points represented over extension fields
  - an explicit matrix report for `π_q` on one chosen basis of `E[n]`,
    together with the congruence checks between matrix trace/determinant mod
    `n` and the global Frobenius invariants
  - explicit Frobenius orbit helpers for both the trivial relative and
    nontrivial absolute finite-field actions
  - fixed points that begin to separate `E(F_p)` from larger finite-field
    point sets represented in the same ambient backend
  - a deterministic visualization layer for those Frobenius objects, including
    traces, characteristic polynomials, local zeta functions, Hasse/type
    reports, extension-count comparisons, characteristic-equation checks,
    torsion/orbit reports, and the twist / isogeny / graph relations
  - a runnable Frobenius milestone tour example that keeps those layers short
    and inspectable from the command line
- the first endomorphism-side arithmetic layer now adds:
  - an educational `QuadraticDiscriminant` value object for integral data such
    as `D = t^2 - 4q`
  - sign and congruence classification for `D < 0` and `D mod 4`
  - a canonical factorization `Δ = v^2 D_K` for negative quadratic-order
    discriminants, with `D_K` fundamental and `v` exposed explicitly
  - an `ImaginaryQuadraticOrder` layer for
    `O_f = Z + f O_K`, with `disc(O_f) = f^2 D_K`
  - relative index computations between included quadratic orders via the
    conductor quotient formula `[O_{f_2} : O_{f_1}] = f_1 / f_2`
  - an `EndomorphismRingCandidateSet` layer that enumerates every candidate
    order between `Z[π]` and `O_K` by listing the divisors `f | v`
  - a Hasse-diagram view of that candidate-order poset, with edge labels
    given by relative indices between immediately adjacent orders
  - an `EndomorphismRingLocalView` layer that records the local conductor gap
    at one prime `\ell` through the valuation `v_\ell(v)`
  - a `VolcanoEndomorphismLevelCandidate` layer that annotates each candidate
    order `O_f` with the arithmetic local level `v_ℓ(f)` without yet claiming
    a certified geometric volcano level for the curve
  - a bridge `EndomorphismVolcanoReport` layer under `isogenies::graphs` that
    compares those arithmetic local candidates with the current graph-theoretic
    volcano heuristic without replacing or overstating that heuristic
  - a tentative `IsogenyEdgeEndomorphismReport` layer that compares the
    source and target local candidate levels of one edge and records only
    possible relations such as horizontal, ascending, or descending
  - automatic node-wise derivation of Frobenius-compatible
    `EndomorphismRingCandidateSet` data directly from stored graph-node curve
    representatives
  - an honest `EndomorphismRingReport` layer that records only the
    Frobenius-compatible sandwich `Z[π] ⊆ End(E) ⊆ O_K` and the remaining
    candidate orders in the ordinary case, while keeping the supersingular
    case separate and explicitly outside the current quadratic-order model
  - fundamental-discriminant detection built on shared exact integer
    squarefreeness helpers in `numerics`
- the first complex-analytic layer adds:
  - upper-half-plane and lattice types
  - truncated Eisenstein sums and analytic invariants
  - the analytic cubic `y² = 4x³ - g₂x - g₃`
  - truncated `℘`, `℘′`, and quasi-periodic `ζ`
  - the torus-to-curve map `z ↦ (℘(z), ℘′(z))`
  - analytic torus torsion and comparison against complex division-polynomial
    `x`-criteria
  - modular `q`-parameters, truncated `q`-expansions for `j`, `E₄`, and `E₆`
  - modular actions, `j`-invariance experiments, and reduction to the
    standard fundamental domain
- the analytic period-recovery and inverse-uniformization layer extends that
  complex-analytic story with:
  period-recovery tooling:
  - Weierstrass cubic-root recovery from `g₂`, `g₃`, including a hybrid
    near-pure-cubic fallback plus Newton polishing
  - cubic-root classification reports and Legendre reduction
  - complex AGM and complete elliptic integral `K`
  - approximate recovery of a period basis, natural `τ`, and canonical `τ`
    modulo `SL₂(ℤ)`
  - validation reports comparing recovered `τ` and lattice invariants against
    the source curve
  - point-level inverse uniformization via an Abel-Jacobi contour integral,
    with roundtrip validation

## Status

This project is still early, but some parts have moved past pure scaffolding.
The focus right now is:

- clear module boundaries
- idiomatic Rust APIs
- mathematical documentation
- tests for basic algebraic behavior
- educational helpers such as textual visualization

The goal is to build a codebase that is pleasant to learn from and easy to
extend correctly.

## What you can study today

- Canonical modular arithmetic through `Fp<P>` and `FpElem<P>`.
- Exact rational arithmetic through `Q`.
- Approximate complex arithmetic through `ComplexApprox`.
- Text-based educational visualization helpers for small prime fields,
  rationals, quotient-polynomial notation, extension-field elements and
  points, complex values, and short Weierstrass curves and points.
- Dense, sparse, and multivariate polynomial representations over fields.
- Univariate evaluation, dense Euclidean division, dense gcd, baseline
  irreducibility classification over `Fp<P>` plus algebraically closed complex
  backends, an exact but partial irreducibility backend for `Q`, and
  first-pass polynomial interpolation via the classical Lagrange formula.
- Real quotient-field modulus checks in `fields::PolynomialModulus` when the
  base field has an irreducibility backend.
- A typed polynomial error surface through `PolynomialError`, shared by
  polynomial arithmetic, evaluation, interpolation, and the corresponding
  explanation helpers.
- Working extension-field arithmetic presented as quotient fields
  `F[x]/(m(x))`, including static tower-friendly designs such as
  `Q(sqrt(2))` and `Q(sqrt(2), i)`.
- Capability traits for extracting square roots, cube roots, and
  characteristic-`p` roots, with the `p`-th-root surface now shared across
  finite-field elements and dense polynomials over finite fields.
- Educational helpers for recurring quadratic extensions such as
  `Fp<P>(sqrt(d))` and `Q(sqrt(d))` through small shared macros that generate
  the corresponding `ExtensionFieldSpec`.
- Operational quotient-value arithmetic through `PolynomialFieldElement<F>`,
  including reduction, quotient-class equality, basic arithmetic, inversion of
  units, and finite-field-oriented explanation helpers.
- Short-Weierstrass curve comparison tools for questions such as:
  - “do these two curves have the same `j`-invariant?”
  - “are they isomorphic over the base field?”
  - “are they only isomorphic after adjoining `sqrt(d)`?”
  - “what is the quadratic twist and how does the scaling map transport points?”
- A first Frobenius surface for short-Weierstrass curves over finite fields,
  including:
  - descriptive metadata for absolute `π_p` and relative `π_q`
  - coefficient-wise Frobenius twists
  - coordinate-wise Frobenius action on affine points
  - Frobenius trace packages recovered from `#E(F_q)` over small finite fields
  - Frobenius discriminants `Δ_π = t^2 - 4q` derived from those trace packages
  - reconstruction helpers between `t` and `#E(F_q)`
  - local zeta functions derived from the Frobenius characteristic polynomial
  - quadratic-twist Frobenius reports derived from the original and twisted traces
  - isogeny Frobenius reports derived from the domain and codomain traces
  - graph-level Frobenius reports derived from isogeny-graph node representatives
  - derived exact counts over `F_{q^n}` from that same trace, both one degree
    at a time and in consecutive prefixes
  - exact Hasse-bound reports derived from the Frobenius trace
  - ordinary versus supersingular classification reports derived from the same trace
  - reports for the relative Frobenius action on exact rational torsion points
  - reports for the absolute Frobenius action on exact torsion points over small
    enumerable extension fields
  - orbit decompositions for relative and absolute Frobenius on rational points
  - small extension-field tests that show how Frobenius fixed points begin to
    distinguish `E(F_p)` from points defined only over larger finite fields
- Division-polynomial tooling for short-Weierstrass curves over
  small enumerable fields, including:
  - explicit base formulas `ψ_0` through `ψ_4`
  - recursive odd/even division-polynomial construction with memoization
  - evaluation at `x` and at affine points
  - rational `x`-candidate recovery, torsion-candidate lifting, exact-order
  torsion filtering, and comparison reports against exhaustive enumeration
  - educational text explanations and a runnable end-to-end example
- Short-Weierstrass function fields and pullback maps, including:
  - the quadratic presentation `F(E) = F(x) ⊕ yF(x)` for one validated short-Weierstrass curve
  - arithmetic, conjugation, norm, inversion, and formal differentiation in that basis
  - pullback maps `phi^* : F(E') -> F(E)` represented by the images of `x'` and `y'`
  - substitution of codomain rational functions and basis elements through those pullbacks
  - contravariant composition of those pullback maps at the function-field level
  - explicit absolute and relative Frobenius isogenies as purely inseparable maps with pullback formulas and differential reports
  - a first Verschiebung surface as a certified function-field-side witness, verified by the relations with Frobenius and `[p]` pullbacks even before full point evaluation is implemented
  - a first certified route to the pullback of `[p]` itself on short-Weierstrass curves, obtained from a Verschiebung certificate rather than from a fully general scalar-multiplication pullback algorithm
- Complex-analytic elliptic-curve tooling, including:
  - validated upper-half-plane points and complex lattices
  - canonical torus coordinates in a fundamental parallelogram
  - truncated Eisenstein sums `G_k`, analytic invariants `g₂`, `g₃`, `Δ`, `j`,
    and the analytic Weierstrass cubic
  - truncated `℘` and `℘′` evaluations with shared approximation metadata
  - the quasi-periodic Weierstrass `ζ` function as a separate analytic surface
  - structured reports for curve membership, the differential equation
    `℘′(z)^2 = 4℘(z)^3 - g₂℘(z) - g₃`, modular invariance, and `q`-expansion
    comparisons
  - torus torsion `E[n] ≅ (1/n)Λ / Λ`, analytic torsion-to-curve mapping, and
    comparison against complex division-polynomial `x`-criteria
  - modular `q = e^{2π i τ}` data, truncated `j(q)`, `E₄(q)`, and `E₆(q)`
  - modular matrices in `SL₂(ℤ)`, their action on `τ`, and reduction to the
    standard fundamental domain
  - high-level analytic aggregate reports such as
    `ComplexAnalyticCurveLabReport` and `UniformizationExperimentReport`
- Analytic period-recovery and inverse-uniformization tooling,
  including:
  - Weierstrass cubic-root recovery from `g₂`, `g₃`
  - cubic-root configuration and recovery reports
  - Legendre reduction, `S₃` orbit inspection, and Legendre-side scale data
  - complex AGM and complete elliptic integrals of the first kind
  - recovery of a full period basis, natural `τ`, and canonical `τ`
  - validation reports based on `j`, recovered lattice invariants, and
    point-level roundtrips `P -> z_P -> (℘(z_P), ℘′(z_P))`
  - reusable lattice-side comparison helpers for approximate equality modulo
    a recovered lattice

## Important caveat

The analytic layer is intentionally educational and approximation-driven.
It is good for studying the structure and for running controlled numerical
experiments, but it is not yet a production-quality complex-analytic engine.

In particular:

- infinite objects are approximated with explicit truncations
- square-box lattice truncations are coordinate-dependent
- convergence quality depends strongly on the chosen `τ`, tolerance, and
  truncation radii
- period recovery and inverse uniformization are intentionally
  educational and heuristic in places; for example, near-equianharmonic
  cubic-root recovery uses a documented hybrid switch rather than a
  production-grade certified solver
- the recovered period basis and recovered torus representatives are
  approximation-driven objects and should be interpreted together with the
  validation reports, not as exact symbolic data

## Examples

The repository now includes concrete examples under:

- [`examples/complex_torus.rs`](./examples/complex_torus.rs)
- [`examples/weierstrass_p.rs`](./examples/weierstrass_p.rs)
- [`examples/fundamental_domain.rs`](./examples/fundamental_domain.rs)
- [`examples/complex_torsion.rs`](./examples/complex_torsion.rs)
- [`examples/root_recovery.rs`](./examples/root_recovery.rs)
- [`examples/legendre-reduction.rs`](./examples/legendre-reduction.rs)
- [`examples/period_recovery.rs`](./examples/period_recovery.rs)
- [`examples/point_roundtrip.rs`](./examples/point_roundtrip.rs)
- [`examples/curve_order.rs`](./examples/curve_order.rs)
- [`examples/frobenius.rs`](./examples/frobenius.rs)
- [`examples/group_structure.rs`](./examples/group_structure.rs)
- [`examples/isomorphism.rs`](./examples/isomorphism.rs)
- [`examples/dual_isogeny.rs`](./examples/dual_isogeny.rs)
- [`examples/division_polynomials.rs`](./examples/division_polynomials.rs)
- [`examples/isogeny_graph.rs`](./examples/isogeny_graph.rs)
- [`examples/velu_isogeny.rs`](./examples/velu_isogeny.rs)
- [`examples/pairing_style_fp12_tower.rs`](./examples/pairing_style_fp12_tower.rs)

Run it with:

```bash
cargo run --example complex_torus
cargo run --example weierstrass_p
cargo run --example fundamental_domain
cargo run --example complex_torsion
cargo run --example root_recovery
cargo run --example legendre-reduction
cargo run --example period_recovery
cargo run --example point_roundtrip
cargo run --example curve_order
cargo run --example frobenius
cargo run --example group_structure
cargo run --example isomorphism
cargo run --example dual_isogeny
cargo run --example division_polynomials
cargo run --example isogeny_graph
cargo run --example velu_isogeny
cargo run --example pairing_style_fp12_tower
```

The pairing-style tower example shows an educational tower

- `Fp`
- `Fp2 = Fp[u] / (u^2 + 1)`
- `Fp6 = Fp2[v] / (v^3 - xi)`
- `Fp12 = Fp6[w] / (w^2 - v)`

Important note:

- the example is intentionally pairing-style, not parameterized for a specific
  production curve such as BLS12-381
- the top tower steps currently use mathematically documented manual
  validation hooks because the crate does not yet expose a generic
  irreducibility backend over arbitrary algebraic-extension bases

The analytic examples are currently:

- `complex_torus`: inspect `τ`, `Λ_τ`, `G₄`, `G₆`, `g₂`, `g₃`, `Δ`, `j`, and
  the analytic versus short-Weierstrass cubic models
- `weierstrass_p`: evaluate `℘(z)` and `℘′(z)` at several torus points and
  check the differential equation numerically
- `fundamental_domain`: reduce a messy `τ` to the standard fundamental domain
  and compare `j` before and after
- `complex_torsion`: study torus torsion, map it to the analytic cubic, and
  compare it against division-polynomial `x`-criteria under the short model
- `root_recovery`: recover the Weierstrass cubic roots from analytic
  invariants and inspect the reconstruction diagnostics
- `legendre-reduction`: compare the six `S₃`-related Legendre parameters and
  inspect what changes versus what stays invariant under root permutation
- `period_recovery`: recover cubic roots, Legendre data, complete elliptic
  integrals, a full period basis, natural `τ`, and canonical `τ`
- `point_roundtrip`: recover a torus representative from a curve point via
  Abel-Jacobi, then validate the roundtrip back through `(℘, ℘′)`

## Example Tour

- Curve order: inspect the order of a short-Weierstrass curve over a small
  prime field through [`examples/curve_order.rs`](./examples/curve_order.rs).
- Frobenius milestone tour: inspect the current finite-field Frobenius layer
  through [`examples/frobenius.rs`](./examples/frobenius.rs).
  The learning goal is to let a reader:
  - start from a trace package `t = q + 1 - #E(F_q)`
  - inspect the derived characteristic polynomial and local zeta function
  - see one concrete pointwise check of the characteristic equation
  - compare one Frobenius-derived extension count against direct enumeration
  - observe a nontrivial absolute-Frobenius orbit over a quadratic extension
  - inspect the twist and isogeny relations through the same trace language
- Group structure: inspect the finite group structure of a small
  short-Weierstrass curve through
  [`examples/group_structure.rs`](./examples/group_structure.rs).
  The example constructs a concrete point, reports the ambient group order,
  prints the order distribution, and summarizes whether the group is cyclic.
- Vélu isogenies: build an educational isogeny lab for small finite elliptic
  curves. The target is to make finite-field isogenies concrete and inspectable
  through explicit kernels, explicit codomain curves, and direct point
  evaluation, starting with small cyclic kernels and Vélu-style formulas over
  small prime fields. The current runnable example lives in
  [`examples/velu_isogeny.rs`](./examples/velu_isogeny.rs).
  The learning goal is to let a reader:
  - start from a short-Weierstrass curve `E / F_p`
  - choose a concrete torsion point and form the finite cyclic subgroup
    generated by it
  - construct the separable isogeny `phi : E -> E'` with that subgroup as
    kernel
  - inspect the resulting codomain curve `E'`
  - evaluate `phi` on explicit points and verify that images land on `E'`
  - observe directly that every kernel point maps to the neutral element and
    that point addition is respected
- Isomorphisms and twists: inspect short-Weierstrass isomorphisms, quadratic twists,
  and base-field versus extension-field comparison through
  [`examples/isomorphism.rs`](./examples/isomorphism.rs).
  The learning goal is to let a reader:
  - compare two curves and see whether they have the same `j`-invariant
  - distinguish “same `j`” from “isomorphic over the current base field”
  - build and explain the explicit scaling map
    `\phi_u(x, y) = (u^2 x, u^3 y)`
  - inspect trivial versus genuinely quadratic twists
  - see a concrete example where a twist becomes isomorphic only after moving
    to `F[x] / (x^2 - d)`
  - inspect the extra automorphisms that appear on the special `j = 1728` and
    `j = 0` loci
- Dual isogenies and composition: compose explicit finite-field isogenies and construct dual
  isogenies by exhaustive search in small examples through
  [`examples/dual_isogeny.rs`](./examples/dual_isogeny.rs).
  The learning goal is to let a reader:
  - compose small explicit isogenies and inspect the multiplicative degree rule
  - treat `[n]` as an explicit educational self-isogeny on `E(F_p)`
  - search exhaustively for a dual of a Vélu isogeny on a tiny finite example
  - verify both dual identities
    `phi_hat ∘ phi = [deg phi]` and `phi ∘ phi_hat = [deg phi]`
    on all enumerated rational points
  - inspect a compact textual summary of the dual construction and the final
    verification checks
- `ℓ`-Isogeny graph exploration: build a small finite-field `ℓ`-isogeny graph explorer for
  short-Weierstrass curves. The learning goal is to let a reader:
  - enumerate rational cyclic kernels of order `ℓ` on a small curve `E / F_p`
  - build the outgoing Vélu isogenies from those kernels
  - deduplicate codomain curves up to base-field isomorphism
  - store the resulting finite graph and inspect its local degree data
  - verify local dual edges using the duality machinery from the dual-isogeny layer
  - inspect educational summaries of connected components, cycles, repeated
    `j`-invariants, and volcano-like heuristic levels
  The current runnable example lives in
  [`examples/isogeny_graph.rs`](./examples/isogeny_graph.rs).
- Division polynomials and rational torsion: implement educational division polynomials for
  short-Weierstrass curves and use them to study rational torsion over small
  finite fields. The learning goal is to let a reader:
  - compute low-degree division polynomials for a concrete curve `E / F_p`
  - evaluate `ψ_n` directly at affine points and compare the odd/even shapes
    `ψ_n ∈ F[x]` and `ψ_n ∈ yF[x]`
  - recover rational `x`-coordinates and rational affine torsion candidates
    from division-polynomial vanishing
  - distinguish raw polynomial candidates from exact-order-`n` torsion points
  - compare polynomial-based torsion detection against explicit point
    enumeration on small finite examples
  - feed recovered torsion generators back into cyclic-kernel and small Vélu
    workflows
  The current runnable example lives in
  [`examples/division_polynomials.rs`](./examples/division_polynomials.rs).
- Complex tori and analytic invariants: build an educational complex-analytic layer for elliptic
  curves over `ℂ`. The learning goal is to let a reader:
  - start from `τ ∈ ℍ` and build the standard lattice `Λ_τ = ℤ + ℤτ`
  - compute truncated Eisenstein sums and derive `g₂`, `g₃`, `Δ`, and `j`
  - pass between the torus `ℂ / Λ` and the analytic cubic
    `y² = 4x³ - g₂x - g₃`
  - evaluate `℘`, `℘′`, and check the differential equation numerically
  - compare analytic `j` computed from lattice sums with `j(q)` from
    `q`-expansions
  - study how `SL₂(ℤ)` acts on `τ` while preserving the underlying complex
    torus, and reduce `τ` to the standard fundamental domain
  - examine torus `n`-torsion analytically and compare it against the complex
    division-polynomial `x`-criterion
  - prepare for later work on period recovery from an analytic cubic
  The current runnable examples live in
  [`examples/complex_torus.rs`](./examples/complex_torus.rs),
  [`examples/weierstrass_p.rs`](./examples/weierstrass_p.rs),
  [`examples/fundamental_domain.rs`](./examples/fundamental_domain.rs), and
  [`examples/complex_torsion.rs`](./examples/complex_torsion.rs).
- Period recovery and inverse uniformization: recover periods and partially invert the uniformization map
  for elliptic curves over `ℂ`. The learning goal is to let a reader:
  - recover the Weierstrass cubic roots from approximate invariants `g₂, g₃`
  - pass from those roots to Legendre form `y² = x(x-1)(x-λ)`
  - evaluate `K(λ)` and `K(1-λ)` through the complex AGM
  - reconstruct a full period basis, the corresponding `τ`, and a canonical
    representative modulo `SL₂(ℤ)`
  - compare the recovered modular data against the source curve using `j` and
    lattice-invariant validation reports
  - partially invert `z ↦ (℘(z), ℘′(z))` at finite points through an
    Abel-Jacobi contour integral and inspect point-level roundtrip residuals
  The current runnable examples live in
  [`examples/root_recovery.rs`](./examples/root_recovery.rs),
  [`examples/legendre-reduction.rs`](./examples/legendre-reduction.rs),
  [`examples/period_recovery.rs`](./examples/period_recovery.rs), and
  [`examples/point_roundtrip.rs`](./examples/point_roundtrip.rs).

## API direction

The library prefers field families that are honest about where their defining
data lives.

Examples:

- prime fields are represented by a compile-time namespace type such as `Fp<17>`
- algebraic extensions are represented by a compile-time field family
  `ExtensionField<S>`, where `S` is an `ExtensionFieldSpec`
- extension-field elements store only their quotient representative value; the
  ambient modulus lives in the specification type, not in each element
- `Field` backends also expose semantic metadata such as
  `IS_ALGEBRAICALLY_CLOSED`, so later APIs can distinguish naturally between
  fields like `Q` and approximate models of `C`

This is intentional: the project tries to keep “what is the field?” separate
from “what is the element?” whenever that makes the math clearer, while still
letting extension fields participate in the same `Field` trait as prime
fields, rationals, and future finite-field towers.

## Visualization philosophy

The project treats explanation helpers as part of the educational API, not as
mere debugging leftovers.

Current visualization helpers focus on deterministic text output, for example:

- operation tables for small prime fields
- step-by-step modular reduction explanations
- exact rational arithmetic explanations
- polynomial formatting in dense, sparse, and multivariate representations
- step-by-step explanations of dense division, dense gcd, polynomial
  evaluation, Lagrange interpolation, dense-polynomial irreducibility, and
  field-modulus irreducibility checks
- readable descriptions of quotient representatives, extension-field arithmetic,
  and modulus suitability
- short-Weierstrass curve and point descriptions, curve-membership checks,
  point-order explanations, scalar-multiplication explanations, and both
  compact and verbose finite-group summaries for small enumerated curves
- short-Weierstrass scaling, isomorphism, quadratic-twist, and
  curve-comparison summaries that make explicit the distinction between
  algebraic-closure isomorphism and base-field isomorphism
- division-polynomial summaries and torsion explanations that report the
  polynomial shape, rational roots, lifted points, exact-order torsion, and
  comparison against exhaustive enumeration
- Frobenius summaries that make explicit the distinction between `π_p` and
  `π_q`, trace packages, characteristic polynomials, local zeta functions,
  Hasse-bound and ordinary/supersingular reports, extension-count comparisons,
  characteristic-equation checks, torsion/orbit reports, and the twist /
  isogeny / graph relations derived from those invariants
- analytic summaries for lattices, Eisenstein truncations,
  analytic invariants, torus-to-curve maps, differential-equation reports,
  modular-action reductions, `q`-expansion comparisons, and analytic torsion
  checks against division polynomials

These helpers are meant to be part of the user-facing learning surface of the
library. They are not just temporary debugging output.

The examples in `examples/` are meant to demonstrate that idea: the objects
are not only constructible, but also inspectable.

If you want step-by-step Mermaid diagrams for some of the core educational
algorithms, see
[`docs/algorithm-diagrams/`](./docs/algorithm-diagrams/README.md).

## Error handling

The repository prefers typed domain errors over ad hoc string failures once an
API starts to stabilize.

Right now that is most visible in `polynomials`, where
[`PolynomialError`](./src/polynomials/error.rs) centralizes recoverable
failures such as:

- division by the zero polynomial
- invalid monic normalization requests on the zero polynomial
- invalid base-field structure for polynomial algorithms
- exact but inconclusive irreducibility attempts in partial backends
- multivariate arity mismatches
- duplicate interpolation abscissas

This is intentional. The educational goal is not only to compute results, but
also to make failure modes mathematically legible.

## Development

Useful commands:

- `cargo fmt`
- `cargo test`
- `cargo clippy --all-targets --all-features`
- `cargo run --example complex_torus`
- `cargo run --example weierstrass_p`
- `cargo run --example fundamental_domain`
- `cargo run --example complex_torsion`
- `cargo run --example curve_order`
- `cargo run --example frobenius`
- `cargo run --example group_structure`
- `cargo run --example isomorphism`
- `cargo run --example dual_isogeny`
- `cargo run --example division_polynomials`
- `cargo run --example isogeny_graph`
- `cargo run --example velu_isogeny`
- `cargo run --example pairing_style_fp12_tower`

## Dependencies

Dependencies are intentionally kept small.

Current numeric dependencies are used narrowly:

- `num-complex` for approximate complex arithmetic
- `num-bigint`, `num-rational`, and `num-traits` for exact arithmetic over `Q`

## Project philosophy

- educational
- correctness before performance
- small explicit APIs
- mathematical assumptions documented in code
- no large algorithmic additions without tests
- exact arithmetic where it improves understanding
- approximation only when clearly labeled
