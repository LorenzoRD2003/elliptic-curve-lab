# AGENTS.md

## Project identity

`elliptic-algorithms-lab` is an educational Rust library for studying and
implementing algebraic, number-theoretic, and cryptographic building blocks.
The repository is intentionally being developed in stages:

- first, clear abstractions
- then, small correct implementations
- then, larger algorithms

This is not a race to maximize features. The main goal is to make the codebase
easy to read, easy to extend, and useful for learning.

## Primary goals

- Build a clean foundation for finite fields, polynomials, elliptic curves, and
  related algorithms.
- Favor explicit mathematical structure over magical APIs.
- Keep the code understandable for someone learning the subject and the Rust
  implementation at the same time.
- Make it easy to inspect, test, and visualize intermediate results.
- Prefer educational output surfaces such as textual explanations, operation
  tables, and polynomial formatting when they help someone understand the math.
- Prefer runnable examples when a design is easier to understand from a small
  end-to-end construction than from API signatures alone.
- Support both finite and infinite base fields when the mathematics naturally
  calls for it, instead of assuming everything is cryptographic or finite from
  the start.
- Expose mathematically meaningful backend properties when they improve later
  APIs, such as whether a field family is algebraically closed.

## Current project posture

- Educational first.
- Correctness before performance.
- Small public APIs before broad feature coverage.
- Step-by-step implementation before optimization.
- Textual explanations and visualizations are welcome when they improve
  understanding.

At the moment, the most mature parts of the repository are `fields` and
`polynomials`, especially:

- `Fp<P>` and `FpElem<P>` for exact prime-field arithmetic
- `Q` for exact rational arithmetic over `BigRational`
- `ComplexApprox` for approximate numerical experiments over `C`
- `ExtensionField<S>` / `ExtensionFieldSpec` as a type-level quotient-field
  design for algebraic extensions and towers over arbitrary base fields,
  including working quotient arithmetic and inversion
- `PolynomialFieldElement<F>` as an autocontained quotient-value layer with
  canonical reduction, quotient-class equality, and basic arithmetic over a
  stored modulus
- `PolynomialModulus<F>::check_field_modulus_requirements()` as the bridge
  from polynomial irreducibility results into field-domain quotient checks
- dense, sparse, and multivariate polynomial representations over fields
- dense Euclidean division and dense gcd over fields
- baseline irreducibility classification over prime fields, plus
  field-theoretic reducibility classification for algebraically closed
  backends such as `ComplexApprox`, plus an exact but partial backend for `Q`
  that certifies some cases and returns an honest inconclusive error otherwise
- univariate evaluation plus baseline interpolation through the classical
  Lagrange formula
- a typed `PolynomialError` surface shared by polynomial-domain APIs and
  explanation helpers
- text-based visualization helpers for prime fields, rationals, polynomials,
  and complex numbers
- runnable educational examples under `examples/`, including extension towers
  that show how the field APIs are meant to be used

## Code style expectations

- Prefer idiomatic, readable Rust over clever or excessively generic code.
- Keep modules small and focused.
- Prefer explicit naming over short cryptic names, especially in educational
  code.
- Public APIs should be conservative and easy to explain.
- Use `Result` for recoverable validation and arithmetic errors.
- Prefer typed domain errors such as `FieldError` and `PolynomialError` over
  raw string errors once a module has more than one meaningful failure mode.
- Add `///` rustdocs to public traits, structs, functions, and any non-obvious
  internal helper that carries important meaning.
- Use `todo!()` only when deferral is intentional and the message explains what
  remains undecided or unimplemented.

## Educational writing rules

- Document mathematical assumptions directly in rustdocs or nearby comments.
- Explain why a representation was chosen when the choice is not obvious.
- Prefer examples and concrete terminology such as `GF(17)` or `F[x]/(m(x))`
  over abstract wording when possible.
- Avoid hiding domain invariants in “smart” helper layers; make them visible in
  types, constructors, or docs.
- If an implementation is approximate, pedagogical, incomplete, or not suitable
  for production cryptography, say so explicitly.

## Architecture conventions

- Keep domain boundaries clear:
  - `fields`: field abstractions and implementations
  - `polynomials`: polynomial representations and later polynomial algorithms
  - `visualization`: educational text-formatting and explanation helpers split
    by mathematical domain
  - `elliptic_curves`: curve models and point representations
  - `algorithms`: reusable algorithmic building blocks
  - `utils`: project-wide helpers that do not belong to a narrower domain
- Re-export only stable, intentional entry points from `lib.rs` and `mod.rs`.
- Prefer lightweight, mathematically honest type-level encodings when they
  remove the need for duplicate runtime context, as in `ExtensionField<S>`.
- Keep error ownership local to the domain:
  - `FieldError` in `fields`
  - `PolynomialError` in `polynomials`
  - avoid duplicating the same failure mode as unrelated strings in several
    files
- When a field family is known at compile time, prefer a namespace type such as
  `Fp<P>`.
- When an algebraic extension can be described statically, prefer a
  specification type plus `ExtensionField<S>` so the extension still
  participates in the main `Field` trait and can itself serve as the base of a
  tower.
- When a higher tower step is mathematically valid but the crate does not yet
  have a generic irreducibility backend for that base field, a documented
  manual validation hook is acceptable in examples and educational extension
  specs. Mark that choice clearly as temporary.
- Do not smuggle ambient field context into element values when a cleaner
  field-family boundary is available.
- Avoid cross-module coupling unless it meaningfully improves clarity.
- Do not add new abstraction layers unless they remove real duplication or
  express a real mathematical boundary.

## Development workflow

Before considering a change complete, run:

- `cargo fmt`
- `cargo test`
- `cargo clippy --all-targets --all-features`

When adding or modifying a runnable example, also run it once if it is cheap
and deterministic.

If a change is intentionally partial, the code should still compile and the
remaining work should be clearly signposted.

## Testing rules

- Do not add large algorithms without tests.
- Prefer deterministic, small examples first.
- For algebraic structures, add property-oriented tests where appropriate:
  - associativity
  - identity laws
  - inverses
  - distributivity
  - compatibility with canonical reduction
- For educational helpers such as formatting or visualization functions, test
  the textual output at the level of important content, not brittle full-file
  snapshots unless the output format is intentionally fixed.
- When a module exposes typed errors, test the error variants directly instead
  of asserting only on formatted messages.
- When polynomial or quotient representations are added, include tests for both
  data invariants and how the chosen storage order is explained to readers.

## Performance rules

- Do not optimize early.
- Prefer the clearest correct implementation first.
- If performance-oriented code is added later, preserve a readable reference
  path when possible.
- Avoid introducing specialized arithmetic tricks, unsafe code, or heavy
  dependencies without a concrete demonstrated need.

## Dependency policy

- No dependency should be added casually.
- A new dependency must have a narrow, justified purpose.
- If a dependency is added, keep the usage small and document why it belongs.
- Prefer standard library facilities unless an external crate materially
  improves correctness, clarity, or educational value.

Current justified numeric dependencies include:

- `num-complex` for approximate complex arithmetic
- `num-bigint` and `num-rational` for exact arithmetic over `Q`
- `num-traits` for numeric identities used by those exact types

## Scope discipline

- This repository is not yet a production cryptography crate.
- Do not present scaffold code as production-safe.
- Do not harden APIs prematurely around features that are not implemented yet.
- Avoid speculative support for serialization, randomness, parallelism, or FFI
  unless the project explicitly moves in that direction.
- Do not assume every algebraic construction should be phrased as a finite
  field. Infinite fields such as `Q` are first-class educational citizens in
  this codebase.

## Final reporting expectations

When summarizing work:

- mention the main files changed
- describe the conceptual change, not just the diff
- note any important simplifications made
- mention remaining risks, TODOs, or intentionally deferred work

The best changes in this repository should feel mathematically honest,
pedagogically useful, and easy for the next contributor to continue.
