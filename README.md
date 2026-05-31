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
  - a scaffold for algebraic extensions through `ExtensionField<F>` and
    `ExtensionFieldDescriptor<F>`
- `polynomials`: early but increasingly structured work on dense, sparse, and
  multivariate representations over fields, with basic arithmetic,
  representation conversions, Euclidean division and gcd in the dense case,
  evaluation, a first interpolation routine, and a shared `PolynomialError`
  surface for recoverable failures.
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
- Univariate evaluation, dense Euclidean division, dense gcd, and first-pass
  polynomial interpolation via the classical Lagrange formula.
- A typed polynomial error surface through `PolynomialError`, shared by
  polynomial arithmetic, evaluation, interpolation, and the corresponding
  explanation helpers.
- The architectural skeleton for extension fields, quotient fields,
  polynomials, and elliptic curves.

## API direction

The library prefers explicit field contexts when runtime metadata matters.

Examples:

- prime fields are represented by a compile-time namespace type such as `Fp<17>`
- extension fields are represented by a runtime field object such as
  `ExtensionField<F>`
- extension-field elements store only their representative value, not the whole
  descriptor of the ambient field

This is intentional: the project tries to keep “what is the field?” separate
from “what is the element?” whenever that makes the math clearer.

## Visualization philosophy

The project treats explanation helpers as part of the educational API, not as
mere debugging leftovers.

Current visualization helpers focus on deterministic text output, for example:

- operation tables for small prime fields
- step-by-step modular reduction explanations
- exact rational arithmetic explanations
- polynomial formatting in dense, sparse, and multivariate representations
- step-by-step explanations of dense division, dense gcd, polynomial
  evaluation, and Lagrange interpolation
- readable descriptions of quotient representatives

These helpers are meant to be part of the user-facing learning surface of the
library. They are not just temporary debugging output.

## Error handling

The repository prefers typed domain errors over ad hoc string failures once an
API starts to stabilize.

Right now that is most visible in `polynomials`, where
[`PolynomialError`](./src/polynomials/error.rs) centralizes recoverable
failures such as:

- division by the zero polynomial
- invalid monic normalization requests on the zero polynomial
- multivariate arity mismatches
- duplicate interpolation abscissas

This is intentional. The educational goal is not only to compute results, but
also to make failure modes mathematically legible.

## Development

Useful commands:

- `cargo fmt`
- `cargo test`
- `cargo clippy --all-targets --all-features`

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
