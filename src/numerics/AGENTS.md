# AGENTS.md for `src/numerics`

## Module mission

The `numerics` module hosts small shared numerical infrastructure that is more
specialized than `utils` but broader than any single mathematical domain.

Its job is to keep approximate policies explicit and reusable across modules
such as `fields` and `elliptic_curves`.

## Design priorities

- Prefer tiny, well-documented value objects.
- Make numerical policy explicit rather than hidden in global constants.
- Keep APIs educational and honest about approximation.
- Avoid turning this module into a grab bag of unrelated helpers.

## Scope guidance

- Good fits include tolerances, truncation settings, normalization choices,
  and similarly small numerical policy types.
- Small exact numerical sequences or coefficient helpers that are shared by
  several mathematical domains may also live here, as long as they are
  documented with their mathematical convention and asymptotic complexity.
- When a classical arithmetic function needs more than one honest evaluation
  strategy, such as divisor-power sums `σ_k(n)`, it is acceptable to keep
  multiple algorithms side by side in `numerics`:
  a literal reference implementation, a better single-input algorithm, and a
  batched `1..=N` algorithm for coefficient-table generation.
- When a shared exact sequence has competing sign conventions, as Bernoulli
  numbers do at `B₁`, document the chosen convention explicitly and normalize
  the implementation to that convention instead of leaving callers to infer it
  from tests.
- For exact arithmetic helpers that feed symbolic coefficient formulas,
  prefer arbitrary-precision integer and rational backends (`BigInt`,
  `BigRational`) over fixed-width integers whenever growth can naturally exceed
  educational toy examples.
- If one consumer needs all values up to a truncation bound, prefer one
  documented batched routine over repeated calls to a single-input helper, and
  say so explicitly in the consuming rustdocs.
- If a type is only meaningful inside one domain and has no shared numerical
  role, keep it local to that domain instead.
- Reexport from consumer modules when that improves ergonomics, but keep the
  canonical definition here when multiple domains depend on it.

## Testing expectations

- Test documented presets directly.
- Test constructor storage behavior directly.
- Keep expected constants explicit in the tests.
- When a tolerance helper implements a mixed absolute/relative comparison,
  test both a near-zero case and a scale-sensitive large-magnitude case.

## Documentation expectations

- State clearly whether a constructor validates inputs or only packages them.
- Say whether a type is exact policy, approximate policy, or experimental
  scaffolding.
- Document the comparison rule explicitly when a helper uses more than a plain
  absolute epsilon.
- For exact-sequence algorithms such as Bernoulli-number generation, document
  the mathematical recurrence or tableau directly and state complexity in
  `Θ(...)` notation, making clear whether the estimate is in arithmetic
  operations, memory, or bit-complexity.
