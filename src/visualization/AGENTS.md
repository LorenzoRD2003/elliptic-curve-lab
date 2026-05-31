# AGENTS.md for `src/visualization`

## Module mission

The `visualization` subtree turns algebraic objects and algorithms into
readable educational text.

Its job is not merely “pretty output”. It should help a reader understand:

1. what mathematical object they are looking at
2. what algorithm or invariant is being explained
3. which parts are exact, approximate, partial, or backend-specific

## Educational posture

- Treat every helper here as teaching material first.
- Prefer short, plain-text explanations over clever formatting.
- Be explicit when a result comes from an approximate backend or from a
  deliberately partial algorithm.
- Avoid pretending an explanation is more general than the underlying code.

## Design priorities

- Clarity before density.
- Stable wording before decorative output.
- Reuse existing formatters and domain helpers instead of duplicating logic.
- Keep explanations close to the mathematical domain they teach.

## Structure rules

- `visualization/fields/` is for field-domain values and field-domain
  explanations.
- `visualization/elliptic_curves/` is for curve equations, points, group-law
  explanations, and small finite curve-group reports.
- `visualization/polynomials/` is for polynomial-domain values and
  polynomial-domain explanations.
- If a helper explains a capability trait such as `SqrtField`, it belongs in
  the matching mathematical subtree rather than in a generic catch-all file.

## Honesty rules

- Say when arithmetic is exact.
- Say when arithmetic is approximate.
- Say when an algorithm is exhaustive but small-scale.
- Say when a curve-group report relies on direct point enumeration or repeated
  addition.
- Say when a backend only handles a subset of mathematically possible cases.
- If a branch choice matters, such as a principal complex square root, say so.

## Formatting guidance

- Prefer deterministic plain text.
- Tables are welcome when they genuinely help, especially for small finite
  fields.
- Explanations should highlight the important intermediate quantities, not
  every possible low-level detail.
- Avoid brittle full-output formatting tricks unless the exact layout is part
  of the teaching goal.

## Testing expectations

- Test important phrases and mathematical content, not giant exact snapshots,
  unless the output format is intentionally fixed.
- For explanatory helpers, test both the “works” case and the “does not exist /
  not supported / no exact answer” case when applicable.
- Keep tests aligned with the actual backend semantics, especially for
  approximate complex arithmetic and exact-only rational helpers.

## Review heuristics

A good change under `src/visualization` should improve at least one of:

- readability
- mathematical honesty
- pedagogical usefulness
- consistency with the underlying domain API

If a visualization helper hides an important caveat from the reader, it is not
done yet.
