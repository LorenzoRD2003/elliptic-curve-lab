use crate::fields::traits::*;
# AGENTS.md for `src/polynomials`

## Module mission

The `polynomials` module is the home for polynomial objects and, later,
polynomial algorithms.

Its role is different from `fields`:

- `fields` models field-like structures and field-specific constructions
- `polynomials` models polynomials as algebraic objects in their own right

This distinction matters. Polynomials should not be treated merely as a helper
for extension fields. They must stay reusable for interpolation, evaluation,
FFT-related work, quotient constructions, and later curve or algorithm layers.

## Educational posture

- This module is educational infrastructure, not a production CAS.
- Prefer code that teaches storage conventions and algebraic meaning clearly.
- Be explicit about whether a method describes raw storage, normalized
  structure, or a true algebraic invariant.
- If something is a scaffold, say so directly in rustdocs.

## Design priorities

- Clarity before cleverness.
- Representation honesty before broad feature coverage.
- Small APIs that are easy to explain.
- Documentation that makes storage order and invariants obvious.
- Reusable polynomial types before advanced algorithms.

## Representation rules

- Dense polynomials should store coefficients in ascending degree order unless a
  different representation is explicitly justified and documented.
- Sparse polynomials should make it obvious which term data is stored and how
  term order is interpreted.
- Multivariate polynomials should document how variable order is chosen and what
  the exponent vector means.
- Do not silently normalize away information unless the API promises that
  canonicalization.
- If trailing-zero trimming is not implemented, say so. If it is implemented,
  say so just as clearly.

## Current scope

At this stage the module is still early, but it is no longer purely
structural. That means:

- representations are welcome
- lightweight metadata helpers are welcome
- educational formatting or explanation helpers are welcome
- small arithmetic and algorithmic building blocks are welcome when they come
  with focused tests
- the `src/visualization/polynomials/` subtree is the right place for
  text-based educational formatting and storage explanations shared across
  polynomial representations
- keep those helpers there instead of re-exporting them from `polynomials::mod.rs`
- keep the `polynomials` root barrel narrow as well:
  - main representation types plus `PolynomialError` are good root candidates
  - term types, irreducibility infrastructure, and shared traits should stay
    under `multivariate`, `sparse`, `irreducibility`, and `traits`
- irreducibility explanations that talk about polynomials as polynomials
  belong in `src/visualization/polynomials/`, not in `fields`

Current practical note:

- `DensePolynomial` is intentionally specialized to coefficient fields via
  `F: Field`
- this is a deliberate scope choice, not a claim that polynomial arithmetic is
  only meaningful over fields in general
- support for more general coefficient rings can be revisited later if the
  project truly needs it
- `DensePolynomial` currently trims trailing zero coefficients so its dense
  storage stays canonical
- `SparsePolynomial` and `MultivariatePolynomial` currently normalize their
  term lists through explicit helper passes that remove zero coefficients and
  combine duplicate monomials
- dense univariate Euclidean division and dense gcd are implemented as the
  current baseline algebraic algorithms over fields
- the shared `UnivariatePolynomial` trait now includes formal derivatives and
  `gcd`
- `SparsePolynomial::gcd` currently delegates through dense conversion to the
  dense Euclidean algorithm instead of maintaining a separate sparse division
  stack
- baseline irreducibility classification is implemented for dense polynomials
  over prime fields, using an intentionally exhaustive educational search,
  and also for algebraically closed backends such as `ComplexApprox` through a
  field-theoretic reducibility conclusion without forced numeric witnesses
- the irreducibility submodule should keep a small public front door:
  - result types and backend capabilities belong in `polynomials::irreducibility`
  - the actual user-facing query should live on `DensePolynomial`, for example
    `polynomial.irreducibility_status()` and `polynomial.is_irreducible()`
  - backend-specific plumbing should stay behind a clearly named internal file
    such as `backend.rs`, not a vague local `traits.rs`
- `Q` now has an exact but partial irreducibility backend:
  - it normalizes to a primitive integer polynomial
  - it searches for small rational roots exactly
  - it uses small-prime Eisenstein checks and irreducibility certificates from
    reductions modulo small primes
  - it returns a typed inconclusive error when those exact criteria do not
    settle the input
- univariate evaluation is implemented for dense and sparse representations
- formal univariate differentiation is implemented for dense and sparse
  representations and exposed through the shared `UnivariatePolynomial` trait
- `p`-th-root extraction over finite fields belongs here when it is a statement
  about polynomials as polynomials, for example the criterion that every
  non-zero term degree must be divisible by the characteristic for a dense
  polynomial to be a `p`-th power in `F[x]`
- Exact integer-root helpers for `IntegerPolynomial` belong here when they are
  small polynomial facts, such as the rational-root theorem over `ℤ[x]`. Curve
  or number-theory consumers should call that helper instead of factoring
  constants and evaluating candidate roots locally. Document their algorithmic
  cost in `Θ(...)` notation, but prefer a readable coarse bound over a highly
  parameterized expression.
- baseline multivariate evaluation is implemented
- Lagrange interpolation is implemented as the first interpolation algorithm
  and should currently live with `DensePolynomial` as a dense-construction
  routine, with other strategies explicitly deferred via TODOs
- `PolynomialError` is the shared failure surface for polynomial-domain APIs
  and should be reused instead of introducing new `&'static str` failures

The current goal is to make dense, sparse, and multivariate polynomial data
structures easy to understand before trying to optimize or generalize them too
far.

## Relationship with `fields`

- Do not move general polynomial concepts into `fields`.
- `fields` may depend on `polynomials` conceptually.
- `polynomials` should remain broader than quotient-field use cases.
- If a type is really a quotient-field or extension-field construction, it
  probably belongs in `fields`.
- If a type is a polynomial representation independent of that use case, it
  belongs here.
- It is acceptable for `polynomials` to depend on the `Field` trait for the
  current stage of development, as long as the docs make that restriction
  explicit and do not overstate generality.
- Field-specific square-root explanations or capability traits belong in
  `fields` / `visualization/fields`, not here, unless a helper is genuinely
  about polynomial roots or polynomial factor behavior.

## Invariants and honesty rules

- Constructors should preserve exactly the invariants they claim to preserve.
- If a constructor performs no normalization, say so.
- If `degree()` reports a stored degree rather than a canonical algebraic
  degree, say so.
- Empty storage is acceptable when the module explicitly treats it as a valid
  zero-polynomial representation.
- Avoid pretending that generic coefficient types support zero-testing,
  canonicalization, or arithmetic unless the type bounds actually guarantee it.
- If a polynomial type currently supports only fields and not general rings,
  say so in the type-level documentation.

## API conventions

- Prefer explicit names such as `coefficients`, `terms`, `arity`, and
  `leading_coefficient`.
- Avoid operator overloading until the arithmetic story is stable enough to be
  pedagogically clear.
- Prefer methods that expose storage meaning before methods that hide it.
- Keep trait bounds minimal unless a method truly needs more structure.
- Do not add large trait hierarchies for polynomial arithmetic yet.
- When an API can fail in a domain-meaningful way, return
  `Result<T, PolynomialError>` instead of an ad hoc string.

## Error conventions

- Keep polynomial-domain recoverable failures in `src/polynomials/error.rs`.
- Prefer specific variants for distinct failure modes, such as:
  - invalid base-field structure for the requested polynomial algorithm
  - zero-polynomial monic normalization
  - division by the zero polynomial
  - multivariate arity mismatch
  - duplicate interpolation abscissas
- Reuse existing variants before adding new ones.
- If visualization helpers simply expose or explain a polynomial algorithm,
  they should usually reuse `PolynomialError` too.
- If an irreducibility backend is intentionally partial, prefer a precise
  `PolynomialError` variant over guessing or silently returning a wrong answer.
- When several consumers need to clear denominators of a `DensePolynomial<Q>`,
  reuse the crate-private rational-normalization helper that produces a
  primitive `Z[x]` representative. Do not duplicate local gcd/lcm/content
  normalization inside downstream algebra modules.
- When an exact integer-root workflow needs a conservative bound from an
  `IntegerPolynomial`, prefer the crate-private Cauchy root bound on that type
  over deriving ad hoc coefficient bounds in downstream modules.

## What not to implement yet

- No heavy generic algebra framework.
- No optimized multiplication strategies yet.
- No FFT-specific polynomial API in the core representations.
- No Gröbner-basis or symbolic-algebra tower unless the project explicitly
  turns in that direction.
- No hidden normalization tricks that make the code harder to explain.

## Testing expectations

Every representation added here should have tests for the properties that make
its storage and semantics understandable:

- storage order
- degree behavior
- empty representation behavior
- leading coefficient behavior
- preservation of constructor input
- normalization behavior, if any
- explanation/formatting output when such helpers are part of the API

When arithmetic is added later, extend tests to cover:

- addition and multiplication on small examples
- formal derivatives, including zero/constant behavior and characteristic-`p`
  cancellation when relevant
- division and gcd where implemented, including monic normalization and zero
  input behavior
- irreducibility classification on small reducible and irreducible examples
- backend-specific irreducibility behavior, including theoretical reducibility
  results without explicit factors when the field metadata justifies them
- exact partial irreducibility behavior for `Q`, including both certified
  answers and honest inconclusive outcomes
- distributivity
- evaluation on small inputs
- interpolation identities where applicable
- agreement between dense and sparse forms when both represent the same object

When evaluation or interpolation code changes, make sure tests cover:

- zero-polynomial behavior
- agreement between equivalent representations
- sample-point reconstruction for interpolation
- rejection paths such as duplicate interpolation abscissas or wrong
  multivariate evaluation arity
- typed error variants when a failure surface is part of the public API

## Documentation expectations

- Public items should explain the mathematical object and the storage model.
- Distinguish “stored degree” from “canonical algebraic degree” when needed.
- If variable order matters, document it.
- If a representation choice is provisional, document the decision point.
- If a caller needs modular exponentiation in `F[x]/(m(x))`, prefer one
  reusable polynomial-side helper such as `pow_mod` over re-implementing
  repeated squaring inside downstream algebra modules.
- Use concrete examples like `[a0, a1, a2]` or `a0 + a1*x + a2*x^2`.
- If a representation participates in the shared polynomial-visualization
  surface, keep that trait small and representation-oriented.
- If a public algorithm can fail, document the main `PolynomialError` cases in
  its rustdocs.

## Review heuristics

A good change to `polynomials` should improve at least one of:

- representation clarity
- mathematical honesty
- reusability
- educational value
- test coverage

If a change makes it harder to answer “what polynomial does this value
represent?”, it is probably moving too fast for the current phase of the
project.
