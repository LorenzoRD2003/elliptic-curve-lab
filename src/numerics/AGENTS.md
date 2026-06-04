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
