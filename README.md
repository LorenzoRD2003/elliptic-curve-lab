# elliptic-algorithms-lab

An educational Rust library for finite fields, polynomials, elliptic curves,
and related algebraic or cryptographic algorithms.

The project is being built in public, in small steps, with an emphasis on:

- mathematical honesty
- readable Rust
- tests for algebraic behavior
- educational formatting and visualization helpers
- architectures that can grow without becoming confusing

Current state after milestone 8:

- milestone 7 is complete enough to study educational division polynomials,
  rational torsion over small finite fields, and their interaction with
  Vélu-style finite-field workflows
- milestone 8 adds the first complex-analytic layer:
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
  - a first period-recovery metadata scaffold for future work

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
  rationals, quotient-polynomial notation, complex values, and short
  Weierstrass curves and points.
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
- Educational helpers for recurring quadratic extensions such as
  `Fp<P>(sqrt(d))` and `Q(sqrt(d))` through small shared macros that generate
  the corresponding `ExtensionFieldSpec`.
- Operational quotient-value arithmetic through `PolynomialFieldElement<F>`,
  including reduction, quotient-class equality, basic arithmetic, inversion of
  units, and prime-field-oriented explanation helpers.
- Short-Weierstrass curve comparison tools for questions such as:
  - “do these two curves have the same `j`-invariant?”
  - “are they isomorphic over the base field?”
  - “are they only isomorphic after adjoining `sqrt(d)`?”
  - “what is the quadratic twist and how does the scaling map transport points?”
- Milestone-7 division-polynomial tooling for short-Weierstrass curves over
  small enumerable fields, including:
  - explicit base formulas `ψ_0` through `ψ_4`
  - recursive odd/even division-polynomial construction with memoization
  - evaluation at `x` and at affine points
  - rational `x`-candidate recovery, torsion-candidate lifting, exact-order
  torsion filtering, and comparison reports against exhaustive enumeration
  - educational text explanations and a runnable end-to-end example
- Milestone-8 complex-analytic elliptic-curve tooling, including:
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
  - a small period-lattice metadata layer for future numerical period recovery

## Important caveat

The analytic milestone is intentionally educational and approximation-driven.
It is good for studying the structure and for running controlled numerical
experiments, but it is not yet a production-quality complex-analytic engine.

In particular:

- infinite objects are approximated with explicit truncations
- square-box lattice truncations are coordinate-dependent
- convergence quality depends strongly on the chosen `τ`, tolerance, and
  truncation radii
- period recovery from `g₂`, `g₃` is not implemented yet; only the metadata
  and validation shape are present

## Examples

The repository now includes concrete examples under:

- [`examples/complex_torus.rs`](./examples/complex_torus.rs)
- [`examples/weierstrass_p.rs`](./examples/weierstrass_p.rs)
- [`examples/fundamental_domain.rs`](./examples/fundamental_domain.rs)
- [`examples/complex_torsion.rs`](./examples/complex_torsion.rs)
- [`examples/curve_order.rs`](./examples/curve_order.rs)
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
cargo run --example curve_order
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

## Milestones

- First milestone: inspect the order of a short-Weierstrass curve over a small
  prime field through [`examples/curve_order.rs`](./examples/curve_order.rs).
- Second milestone: inspect the finite group structure of a small
  short-Weierstrass curve through
  [`examples/group_structure.rs`](./examples/group_structure.rs).
  The example constructs a concrete point, reports the ambient group order,
  prints the order distribution, and summarizes whether the group is cyclic.
- Third milestone: build an educational isogeny lab for small finite elliptic
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
- Fourth milestone: inspect short-Weierstrass isomorphisms, quadratic twists,
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
- Fifth milestone: compose explicit finite-field isogenies and construct dual
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
- Sixth milestone: build a small finite-field `ℓ`-isogeny graph explorer for
  short-Weierstrass curves. The learning goal is to let a reader:
  - enumerate rational cyclic kernels of order `ℓ` on a small curve `E / F_p`
  - build the outgoing Vélu isogenies from those kernels
  - deduplicate codomain curves up to base-field isomorphism
  - store the resulting finite graph and inspect its local degree data
  - verify local dual edges using the milestone-5 duality machinery
  - inspect educational summaries of connected components, cycles, repeated
    `j`-invariants, and volcano-like heuristic levels
  The current runnable example lives in
  [`examples/isogeny_graph.rs`](./examples/isogeny_graph.rs).
- Seventh milestone: implement educational division polynomials for
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
- Eighth milestone: build an educational complex-analytic layer for elliptic
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
  - prepare for a future milestone on period recovery from an analytic cubic
  The current runnable examples live in
  [`examples/complex_torus.rs`](./examples/complex_torus.rs),
  [`examples/weierstrass_p.rs`](./examples/weierstrass_p.rs),
  [`examples/fundamental_domain.rs`](./examples/fundamental_domain.rs), and
  [`examples/complex_torsion.rs`](./examples/complex_torsion.rs).

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
- milestone-8 analytic summaries for lattices, Eisenstein truncations,
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
