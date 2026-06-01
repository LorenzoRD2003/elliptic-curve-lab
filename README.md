# elliptic-algorithms-lab

An educational Rust library for finite fields, polynomials, elliptic curves,
and related algebraic or cryptographic algorithms.

The project is being built in public, in small steps, with an emphasis on:

- mathematical honesty
- readable Rust
- tests for algebraic behavior
- educational formatting and visualization helpers
- architectures that can grow without becoming confusing

## Status

This project is still early, but some parts have moved past pure scaffolding.
The focus right now is:

- clear module boundaries
- idiomatic Rust APIs
- mathematical documentation
- tests for basic algebraic behavior
- educational helpers such as textual visualization

The current center of gravity is `src/fields/` plus the core pieces of
`src/polynomials/`.

The goal is not to ship every algorithm immediately. The goal is to build a
codebase that is pleasant to learn from and easy to extend correctly.

## Current areas

- `fields`: the most developed module so far. It includes:
  - exact prime fields through `Fp<P>` and `FpElem<P>`
  - exact rationals through `Q`
  - approximate complex arithmetic through `ComplexApprox`
  - statically specified algebraic extensions through
    `ExtensionField<S>` and `ExtensionFieldSpec`, with quotient reduction,
    multiplication, inversion, and tower-friendly composition
- `polynomials`: early but increasingly structured work on dense, sparse, and
  multivariate representations over fields, with basic arithmetic,
  representation conversions, Euclidean division and gcd in the dense case,
  baseline irreducibility checks over prime fields plus algebraically closed
  backends such as `ComplexApprox`, a first exact but partial backend over
  `Q`, evaluation, a first interpolation routine, and a shared
  `PolynomialError` surface for recoverable failures.
- `visualization`: shared educational text-formatting and explanation helpers,
  now covering `fields`, `polynomials`, and the first elliptic-curve group
  summaries and explanations under `visualization/elliptic_curves/`.
- `elliptic_curves`: early structural scaffolding.
- `algorithms`: placeholder modules for future reusable algorithms.

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
- Operational quotient-value arithmetic through `PolynomialFieldElement<F>`,
  including reduction, quotient-class equality, basic arithmetic, inversion of
  units, and prime-field-oriented explanation helpers.
- The architectural skeleton for quotient fields and elliptic curves.

## Examples

The repository now includes concrete examples under
[`examples/curve_order.rs`](./examples/curve_order.rs),
[`examples/group_structure.rs`](./examples/group_structure.rs),
[`examples/velu_isogeny.rs`](./examples/velu_isogeny.rs),
and a larger extension-field example under
[`examples/pairing_style_fp12_tower.rs`](./examples/pairing_style_fp12_tower.rs).

Run it with:

```bash
cargo run --example curve_order
cargo run --example group_structure
cargo run --example velu_isogeny
cargo run --example pairing_style_fp12_tower
```

```rust
use elliptic_algorithms_lab::{EnumerableCurveModel, Field, Fp, ShortWeierstrassCurve};

type F = Fp<101>;

let e = ShortWeierstrassCurve::<F>::new(F::from_i64(2), F::from_i64(3))?;
let points = e.points();
let order = e.order();

println!("#E(F_101) = {}", order);
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
- that is a teaching choice, not a claim of production-ready pairing-field
  infrastructure

## Milestones

- First milestone: inspect the order of a short-Weierstrass curve over a small
  prime field through [`examples/curve_order.rs`](./examples/curve_order.rs).
- Second milestone: inspect the finite group structure of a small
  short-Weierstrass curve through
  [`examples/group_structure.rs`](./examples/group_structure.rs),
  where the library constructs a concrete point, reports the ambient group
  order, prints the order distribution, and summarizes whether the group is
  cyclic:
- Third milestone: build an educational isogeny lab for small finite elliptic
  curves. The target is to make finite-field isogenies concrete and inspectable
  through explicit kernels, explicit codomain curves, and direct point
  evaluation, starting with small cyclic kernels and Vélu-style formulas over
  small prime fields. The current runnable example lives in
  [`examples/velu_isogeny.rs`](./examples/velu_isogeny.rs).
  In mathematical terms, the learning goal is to let a reader:
  - start from a short-Weierstrass curve `E / F_p`
  - choose a concrete torsion point and form the finite cyclic subgroup
    generated by it
  - construct the separable isogeny `phi : E -> E'` with that subgroup as
    kernel
  - inspect the resulting codomain curve `E'`
  - evaluate `phi` on explicit points and verify that images land on `E'`
  - observe directly that every kernel point maps to the neutral element and
    that point addition is respected

  The current milestone example uses the curve `E / F_101` given by
  `y^2 = x^3 + 2x + 3`, the order-`3` point `P = (35, 15)`, and the
  corresponding cyclic kernel `<P>`. Running
  `cargo run --example velu_isogeny` prints the kernel, the computed codomain,
  and step-by-step explanations of the Vélu codomain formulas and the
  point-evaluation formulas.

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

These helpers are meant to be part of the user-facing learning surface of the
library. They are not just temporary debugging output.

The examples in `examples/` are meant to demonstrate that idea: the objects
are not only constructible, but also inspectable.

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
- `cargo run --example curve_order`
- `cargo run --example group_structure`
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
