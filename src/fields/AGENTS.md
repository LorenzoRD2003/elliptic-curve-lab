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
- Separate field context from element values unless there is a strong reason
  not to.

## Educational posture

- Treat this module as educational infrastructure, not as production crypto.
- Say when something is exact, approximate, scaffold-only, or intentionally
  deferred.
- Favor functions that expose structure and intermediate meaning over hiding
  everything behind overloaded traits or macros.
- When a helper is meant for pedagogy, visualization, or explanation, that is
  a feature, not noise.

## Core invariants

- Prime-field elements should store canonical representatives whenever the type
  promises that behavior.
- Rational values should stay exact and normalized through `BigRational`
  instead of introducing approximate shortcuts.
- Constructors that return `Result` should enforce the invariants they claim to
  enforce.
- If a function assumes a modulus or extension specification was already
  validated, document that assumption clearly.
- Keep the distinction clear between:
  - structural modulus construction through `PolynomialModulus::new`
  - stronger quotient-field validation through
    `PolynomialModulus::check_field_modulus_requirements`
- If a field backend exposes semantic mathematical metadata such as algebraic
  closedness, document whether that statement refers to an exact field or to an
  approximate numerical model of one.
- Non-zero mathematical requirements should be expressed in the type system when
  it improves clarity, such as `NonZeroU32` for extension degrees.
- Distinct field specifications or incompatible parameterizations must not be
  mixed silently.
- Extension-field elements should not silently carry ambient runtime metadata.
  The current direction is:
  - `ExtensionFieldSpec` for static metadata and validation
  - `ExtensionField<S>` for the field family itself
  - `ExtensionFieldElement<S>` for quotient representatives only
  - when those quotient representatives are stored canonically, it is
    acceptable to implement structural traits such as `Hash` directly from the
    trimmed coefficient vector so higher layers can build small hashed lookup
    tables without reintroducing runtime field metadata

## Trait conventions

- `Field` should stay focused on the smallest useful algebraic interface:
  identities, arithmetic, inversion, equality, and simple embedding helpers.
- When a mathematically natural field family depends on runtime ambient data
  rather than on type-level data alone, prefer a separate trait such as
  `AmbientField` over forcing that family into the static `Field` interface.
- Semantic field-family metadata is welcome in `Field` when it captures a real
  mathematical property that later APIs can build on. The current examples are
  `IS_ALGEBRAICALLY_CLOSED` and field characteristic.
- `FiniteField` should cover field metadata and structural checks, not every
  possible algorithm over finite fields.
- Capability traits such as `SqrtField` are encouraged when an operation is
  real, useful, and only honestly implementable for some backends.
- The same applies to characteristic-`p` capabilities such as
  `PthRootExtraction`:
  - finite-field elements admit a canonical `p`-th root because Frobenius is
    invertible
  - broader algebraic objects such as polynomials may admit a `p`-th root only
    on a smaller locus
  - keep those distinctions explicit instead of pretending every
    characteristic-`p` object has the same extraction story
  - do not reuse `PthRootExtraction` for “invert the pullback of Frobenius on a
    twisted ambient object” when the mathematical operation is really a
    coordinate-substitution inverse rather than an honest `p`-th root in the
    same algebraic structure
- `EnumerableFiniteField` is also an acceptable capability trait when the
  backend can honestly enumerate all elements and that enumeration is still
  educationally reasonable.
- `SqrtField` should remain small and honest:
  - it should promise only square-root discovery, not a full quadratic-solving
    framework
  - its docs should state when an implementation is exhaustive, exact-only,
    approximate, branch-sensitive, or intentionally partial
  - for small finite extension backends, a documented brute-force square-root
    search over an honestly enumerable field is acceptable as an educational
    first implementation
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
  - `NonInvertibleElement`
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
- Keep extension specifications lightweight and descriptive.
- It is fine for quotient values such as `PolynomialFieldElement<F>` to be
  operational and self-contained without also becoming the primary field-family
  backend of the crate.
- The same is now true for `RationalFunction<F>`:
  - keep it as an autocontained canonical value type first
  - use `DensePolynomial<F>` as the current honest backend
  - prefer eager gcd reduction plus denominator-monic normalization over
    storing arbitrary presentations
  - value-intrinsic queries such as “is this rational function constant?”
    belong on `RationalFunction<F>` itself rather than in downstream
    consumers
- When the rational-function layer needs field-family integration, prefer a
  separate zero-sized `RationalFunctionField<F>` whose `Elem` is
  `RationalFunction<F>` instead of collapsing family metadata and stored
  value into one type
- Prefer explicit element constructors over implicit conversions when invariants
  matter.
- Avoid internal representations that are hard to explain unless there is a
  strong payoff.
- When an algorithm has meaningful backend-specific interpretation, keep that
  honesty visible in nearby docs. Current examples:
  - `Fp<P>` square roots via Tonelli-Shanks for odd primes
  - `Q` square roots only when the rational is already a square in `Q`
  - `ComplexApprox` square roots as approximate principal-branch values
  - shared numerical tolerance policy should live in sibling infrastructure
    such as `src/numerics/` rather than in `elliptic_curves`

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
- No claim that every field backend supports square roots just because
  `SqrtField` exists.
- No claim that every finite field backend should enumerate all elements just
  because `EnumerableFiniteField` exists.
- No premature general square-root framework for arbitrary extension fields
  until the crate really needs it and can explain it clearly.
  A small finite-extension implementation may still use exhaustive
  enumeration when the docs say clearly that the algorithm is educational and
  only appropriate for tiny fields.

## Visualization and debugging guidance

- Educational visualization belongs here when it helps explain the math.
- Prefer text-based, deterministic formatting and explanation helpers first.
- Functions that build tables, explain reductions, or format elements clearly
  are encouraged.
- When finite-field family metadata needs a human-readable mathematical label,
  prefer keeping that formatting on `FiniteFieldDescriptor` itself, so
  downstream modules reuse one shared notation for `F_p` and `F_{p^n}` rather
  than inventing local string conventions.
- Runnable examples that show how a field family is assembled and used are also
  encouraged when they make the abstraction easier to learn.
- Keep those explanatory helpers under `src/visualization/fields/`; do not
  re-export them from `fields::mod.rs` just to shorten imports.
- For infinite fields such as `Q`, prefer explanations of exact arithmetic,
  canonical forms, inverses, and quotient notation over impossible or
  misleading “full tables”.
- Square-root explanations should state clearly whether the backend is exact,
  approximate, exhaustive, or algorithmic.
- Element-enumeration helpers should state clearly that they are intended only
  for small finite educational backends.
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

For algebraic extensions, prefer designs that keep the quotient presentation at
the type level when that enables towers and lets the backend implement the main
`Field` trait directly.

When the extension-field layer grows beyond one concern, prefer splitting it as
its own `extension_field/` module tree instead of keeping macros, quotient
types, trait impls, finite-field capabilities, and tests inside one long file.
The current preferred separation is:

- `mod.rs` for the public surface and shared internal aliases
- `macros.rs` for educational quadratic-extension helpers
- `spec.rs` for `ExtensionFieldSpec`
- `field.rs` for the zero-sized field family and quotient arithmetic helpers
- `element.rs` for stored quotient representatives
- `traits.rs` for `Field`/`FiniteField`/capability impls
- `tests.rs` for module-local regression coverage

The same separation heuristic applies to `rational_function_field/` now that it
contains both:

- `value.rs` for the canonical stored rational-function value
- `field.rs` for the zero-sized field family and `Field` impl
- `tests.rs` for module-local regression coverage

For extension towers used as examples:

- it is acceptable to use mathematically documented manual validation hooks for
  upper tower steps when generic irreducibility support over the intermediate
  base field does not exist yet
- say so explicitly in code comments and user-facing docs
- prefer examples that teach the tower shape clearly over examples that pretend
  to be production pairing-field parameter sets
- for recurring quadratic examples over `Fp<P>`, prefer the shared
  `define_fp_quadratic_extension!` helper over rewriting identical
  `ExtensionFieldSpec` boilerplate in each test module
- for recurring `Q(sqrt(d))` examples, prefer the shared
  `define_q_quadratic_extension!` helper over rewriting identical rational
  quadratic-extension specs in each test module

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
- if a square-root capability is implemented:
  - roots square back to the input
  - non-residues or unsupported cases are rejected honestly
  - paired roots are additive inverses when that is the promised behavior

For exact fields such as `Q`, also test:

- canonical normalization behavior
- exact inverse and division behavior
- integer embeddings
- exact `pow` / `square` / `cube` behavior on small examples

Use small finite examples whenever possible. `Fp<17>`-style tests are
excellent for teaching and for catching regressions.

## Documentation expectations

- Public items should explain what they model, not just what they return.
- If a method assumes validated structure, say so.
- If equality is approximate, say so.
- If arithmetic is exact, say so too.
- If approximate arithmetic uses a shared tolerance object, document both the
  default preset and the explicit override path.
- If an approximate backend exposes comparison reports, store the actual
  compared element values and the tolerance policy explicitly in that report.
- If a function is intentionally educational rather than optimal, say so.
- If an algorithm is complete only on a subset of inputs or only exact for a
  subset of the mathematical field, say so directly where the API lives.
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
