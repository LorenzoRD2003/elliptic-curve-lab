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

## First milestone

The goal of the first milestone is to make the following kind of example
trivial for a user to write and understand:

```rust
type F = Fp<101>;

let e = ShortWeierstrassCurve::<F>::new(F::from(2), F::from(3))?;
let points = e.points();
let order = e.order_naive();

println!("#E(F_101) = {}", order);
```

That example captures the intended learning experience well:

- choose a concrete field
- define a curve over that field
- enumerate or inspect its points
- compute a naive order in a way that is easy to read

The milestone is not just "have the APIs exist". It is to make the path from
the mathematics to runnable code feel direct, unsurprising, and educational.

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
  with `visualization/fields/` and `visualization/polynomials/` as sibling
  branches.
- `elliptic_curves`: early structural scaffolding.
- `algorithms`: placeholder modules for future reusable algorithms.

## What you can study today

- Canonical modular arithmetic through `Fp<P>` and `FpElem<P>`.
- Exact rational arithmetic through `Q`.
- Approximate complex arithmetic through `ComplexApprox`.
- Text-based educational visualization helpers for small prime fields,
  rationals, quotient-polynomial notation, and complex values.
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

The repository now includes a concrete example under
[`examples/pairing_style_fp12_tower.rs`](./examples/pairing_style_fp12_tower.rs).

Run it with:

```bash
cargo run --example pairing_style_fp12_tower
```

That example shows an educational tower

- `Fp`
- `Fp2 = Fp[u] / (u^2 + 1)`
- `Fp6 = Fp2[v] / (v^3 - xi)`
- `Fp12 = Fp6[w] / (w^2 - v)`

with readable textual output for:

- each extension presentation
- the tower generators `u`, `v`, and `w`
- the defining quotient relations after reduction
- a sample multiplication trace inside `Fp12`

Important note:

- the example is intentionally pairing-style, not parameterized for a specific
  production curve such as BLS12-381
- the top tower steps currently use mathematically documented manual
  validation hooks because the crate does not yet expose a generic
  irreducibility backend over arbitrary algebraic-extension bases
- that is a teaching choice, not a claim of production-ready pairing-field
  infrastructure

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

These helpers are meant to be part of the user-facing learning surface of the
library. They are not just temporary debugging output.

The extension-field example in `examples/` is meant to demonstrate that idea:
the tower is not only constructible, but also inspectable.

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
