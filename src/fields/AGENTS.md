# AGENTS.md for `src/fields`

## Module mission

The `fields` module is the mathematical foundation of the library. It should
provide clean, well-documented abstractions for field-like structures and a
small set of concrete implementations that are easy to study.

This module has two jobs:

1. express the algebra clearly
2. teach the reader what the implementation is doing

Code here should be understandable to someone learning finite fields, modular
arithmetic, and Rust at the same time.

That mission now includes both exact finite-field examples and exact infinite
base fields such as `Q`, plus one clearly labeled approximate field-like model
for numerical intuition.

## Design priorities

- Clarity before cleverness.
- Correct invariants before broad API surface.
- Canonical representations where practical.
- Rich documentation for public APIs.
- Small, reviewable steps.
- Separate field context from element values unless there is a strong reason not
  to.

## Educational posture

- Treat this module as educational infrastructure, not as production crypto.
- Say when something is exact, approximate, scaffold-only, or intentionally
  deferred.
- Favor functions that expose structure and intermediate meaning over hiding
  everything behind overloaded traits or macros.
- When a helper is meant for pedagogy, visualization, or explanation, that is a
  feature, not noise.

## Core invariants

- Prime-field elements should store canonical representatives whenever the type
  promises that behavior.
- Rational values should stay exact and normalized through `BigRational`
  instead of introducing approximate shortcuts.
- Constructors that return `Result` should enforce the invariants they claim to
  enforce.
- If a function assumes a field descriptor or modulus was already validated,
  document that assumption clearly.
- Keep the distinction clear between:
  - structural modulus construction through `PolynomialModulus::new`
  - stronger quotient-field validation through
    `PolynomialModulus::check_field_modulus_requirements`
- If a field backend exposes semantic mathematical metadata such as algebraic
  closedness, document whether that statement refers to an exact field or to an
  approximate numerical model of one.
- Non-zero mathematical requirements should be expressed in the type system when
  it improves clarity, such as `NonZeroU32` for extension degrees.
- Distinct field descriptors or incompatible parameterizations must not be mixed
  silently.
- Extension-field elements should not silently carry ambient runtime metadata
  unless the design has been intentionally revisited. The current direction is:
  - `ExtensionFieldDescriptor<F>` for metadata
  - `ExtensionField<F>` for the ambient runtime field object
  - `ExtensionFieldElement<F>` for value representatives only

## Trait conventions

- `Field` should stay focused on the smallest useful algebraic interface:
  identities, arithmetic, inversion, equality, and simple embedding helpers.
- Semantic field-family metadata is welcome in `Field` when it captures a real
  mathematical property that later APIs can build on. The current example is
  `IS_ALGEBRAICALLY_CLOSED`.
- `FiniteField` should cover field metadata and structural checks, not every
  possible algorithm over finite fields.
- Not every `Field` implementor should also be a `FiniteField`; keep that
  distinction meaningful.
- Do not overload core traits with serialization, randomness, FFT hooks, or
  unrelated conveniences unless the project explicitly needs them.
- Prefer one trait with a clear conceptual boundary over many tiny traits that
  fragment the API without improving understanding.

## Error conventions

- Use `FieldError` for domain-level failures.
- Keep `FieldError` local to field-domain failures and do not reuse it for
  polynomial-only concerns now that `polynomials` has its own
  `PolynomialError`.
- Prefer specific variants such as:
  - `DivisionByZero`
  - `InvalidModulus`
  - `InvalidPolynomialModulus`
  - `NonIrreduciblePolynomial`
  - `ElementOutOfRange`
  - `IncompatibleFieldParameters`
- Add a new error variant only if it reveals a distinct mathematical or API
  failure mode.
- Do not hide meaningful failure reasons behind generic strings unless the case
  is genuinely temporary scaffolding.

## Representation rules

- For prime fields, prefer small explicit types with transparent invariants.
- For rationals, prefer exact arbitrary-precision arithmetic over hand-rolled
  approximations.
- For extension and quotient fields, preserve room to evolve the arithmetic
  representation while keeping the API conceptually honest.
- Keep descriptors lightweight and descriptive.
- Prefer explicit element constructors over implicit conversions when invariants
  matter.
- Avoid internal representations that are hard to explain unless there is a
  strong payoff.

## What not to implement yet

- No optimized Montgomery, Barrett, or similar reduction machinery.
- No advanced irreducibility testing framework.
- It is fine, however, to integrate the currently available polynomial
  irreducibility backends into `fields` when the goal is to validate quotient
  moduli honestly.
- No giant polynomial arithmetic subsystem unless it arrives with focused tests.
- No FFT- or pairing-specific field hooks in the base abstractions.
- No production-style trait explosion for every conceivable algebraic nuance.
- No unsafe code for field arithmetic at this stage.

## Visualization and debugging guidance

- Educational visualization belongs here when it helps explain the math.
- Prefer text-based, deterministic formatting and explanation helpers first.
- Functions that build tables, explain reductions, or format elements clearly
  are encouraged.
- For infinite fields such as `Q`, prefer explanations of exact arithmetic,
  canonical forms, inverses, and quotient notation over impossible or
  misleading “full tables”.
- Polynomial helpers that explain coefficient order, quotient notation, or
  modulus role are encouraged and fit the educational mission well.
- Avoid pulling in plotting or graphical dependencies unless there is a clear,
  concrete educational gain.

## Adding a new field implementation

When adding a new field family:

1. create one focused file under `src/fields/`
2. define the representation and document it
3. explain which invariants are enforced by constructors
4. wire the module through `src/fields/mod.rs` only if the API is coherent
5. add tests for arithmetic behavior and invariants
6. add educational rustdocs, especially if the field is approximate or unusual

If the field also deserves explanatory helpers, add a matching file under
`src/visualization/fields/` and wire it through the public reexports only when
the API is coherent and stable enough to teach from.

## Testing expectations

Any new algebraic implementation or algorithm under `fields` should be backed
by tests for the properties that apply:

- associativity
- additive identity
- multiplicative identity
- additive inverses
- multiplicative inverses when defined
- distributivity
- canonical reduction behavior
- rejection of invalid construction inputs
- storage-order expectations for polynomial coefficients
- educational formatting or visualization output when a helper is meant to
  explain structure rather than just compute it

For exact fields such as `Q`, also test:

- canonical normalization behavior
- exact inverse and division behavior
- integer embeddings
- exact `pow` / `square` / `cube` behavior on small examples

Use small finite examples whenever possible. `Fp<17>`-style tests are excellent
for teaching and for catching regressions.

## Documentation expectations

- Public items should explain what they model, not just what they return.
- If a method assumes validated structure, say so.
- If equality is approximate, say so.
- If arithmetic is exact, say so too.
- If a function is intentionally educational rather than optimal, say so.
- If a future design decision is unresolved, record the decision point in a
  short and concrete `todo!()` or rustdoc note.

## Review heuristics

A good change to `fields` should make at least one of these better:

- mathematical honesty
- readability
- invariant safety
- educational value
- test coverage

If a change makes the code harder to explain than before, it is probably moving
too fast for this stage of the project.
