# elliptic-algorithms-lab

An educational Rust library for finite fields, polynomials, elliptic curves,
and adjacent algebraic tooling.

The project is intentionally not trying to be a production cryptography crate.
The goal is to build small, readable mathematical components that are honest
about their scope, easy to test, and pleasant to learn from.

## What the crate is trying to do

- build clear foundations for field arithmetic, polynomial arithmetic, and
  elliptic curves
- expose mathematically meaningful APIs instead of hiding everything behind
  generic wrappers
- keep explanations, reports, and examples as first-class educational surfaces
- grow in small steps, with correctness and readability ahead of coverage or
  performance

## What is already usable

Today the repo is most useful as a lab for:

- prime fields, rationals, algebraic extensions, and rational-function fields
- dense, sparse, and multivariate polynomials
- short-Weierstrass curves over small fields
- unified curve-side group-order selection across exhaustive,
  quadratic-character, and prime-field Mestre routes
- short-Weierstrass group-order parity from the rational `2`-torsion gcd
  criterion
- point-order recovery from one known annihilating multiple
- naive Hasse-interval search for an annihilating multiple `[M]P = O`
- unified curve-side point-order selection across exhaustive, known-multiple,
  and naive Hasse-interval routes
- unified curve-side group-exponent recovery across exhaustive and sampled
  point-order accumulation routes
- separate verification of one sampled exponent lower bound against the
  Hasse interval coming from a chosen group-order route
- torsion, division polynomials, and small explicit isogeny workflows
- Frobenius data over finite fields, including traces, characteristic
  polynomials, character-sum counts, Hasse checks, Hasse intervals, extension
  counts, and related reports
- a first endomorphism-side layer derived from Frobenius discriminants
- a substantial complex-analytic layer around lattices, `℘`, modular data,
  period recovery, and inverse uniformization
- deterministic text visualization helpers that explain what the library is
  computing

## Project style

The repo prefers:

- educational honesty over polished marketing claims
- explicit types and narrow capability traits
- small modules with local error ownership
- exact arithmetic where it is natural, and explicit approximation reports when
  it is not

If something is approximate, tiny-field-only, exhaustive, heuristic, or still
scaffold-level, the library tries to say so directly.

## Where to start

If you want to explore the crate from the command line, these examples are good
entry points:

- `cargo run --example curve_order`
- `cargo run --example frobenius`
- `cargo run --example division_polynomials`
- `cargo run --example velu_isogeny`
- `cargo run --example isogeny_graph`
- `cargo run --example complex_torus`
- `cargo run --example period_recovery`
- `cargo run --example point_roundtrip`

These cover the main finite-field, isogeny, Frobenius, and analytic threads of
the project without forcing a full tour of every feature.

## Caveats

- The finite-field elliptic-curve side is strongest on small enumerable
  examples and educational reports.
- The analytic side is approximation-driven and intended for controlled
  experiments, not production numerics.
- Some advanced surfaces are intentionally partial: the crate prefers a narrow
  honest implementation over a broad misleading one.

## Documentation philosophy

The README is intentionally brief. The detailed mathematical story is meant to
live in:

- module docs and rustdocs
- typed reports and visualization helpers
- runnable examples under `examples/`
- the repository guidance in `AGENTS.md`

If you want to extend the library, start there rather than treating this file
as a full feature inventory.
