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
  explanations, small finite curve-group reports, and milestone-7
  division-polynomial / rational-torsion summaries.
- `visualization/isogenies/` is for educational summaries of kernels,
  codomains, point-evaluation formulas for explicit isogeny constructions, and
  milestone-level summaries of composition, scalar multiplication, dual
  isogenies, exhaustive dual verification on tiny curves, and milestone-6
  graph summaries and adjacency explanations.
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
- Say when a division-polynomial explanation relies on exhaustive root scans,
  point lifting, or exact-order filtering on small finite fields.
- Say when an isogeny explanation is tied to a specific normalized formula such
  as the current short-Weierstrass Vélu construction.
- Say when the explanation is only valid for explicit finite kernels rather
  than arbitrary subgroup schemes.
- Say when a backend only handles a subset of mathematically possible cases.
- If a branch choice matters, such as a principal complex square root, say so.

## Formatting guidance

- Prefer deterministic plain text.
- Tables are welcome when they genuinely help, especially for small finite
  fields.
- Explanations should highlight the important intermediate quantities, not
  every possible low-level detail.
- For isogenies, prefer showing the kernel points, the codomain formulas, and
  a few key translation-sum terms over dumping large algebraic expressions
  without guidance.
- For graph visualizations, prefer one compact structural summary plus explicit
  node/edge listings and adjacency lists. Say directly that nodes store
  representatives and edges may carry transport witnesses onto those stored
  representatives.
- For milestone-7 division-polynomial explanations, prefer showing:
  - the curve and the index `n`
  - the shape of `ψ_n`
  - the polynomial obtained
  - rational roots and lifted points
  - torsion / exact-order filtering results
  - comparison against exhaustive enumeration
- If a graph summary includes a volcano-like layering, keep it explicitly
  heuristic and explain how its root was chosen. A deterministic weak-component
  root plus role counts is acceptable; do not present it as arithmetic proof of
  a true isogeny volcano.
- Prefer names in visualization that preserve that boundary in plain sight, for
  example `VolcanoLike...` rather than a name that sounds like certified
  structure.
- Do not let volcano-like presentation drift from visual intuition into
  mathematical-sounding certification. If wording starts to imply theorem-level
  structure, either add the missing mathematical justification or soften the
  explanation back to a clearly pedagogical heuristic.
- If a separate helper explains one inferred layering, it should reuse the
  already computed levels and roles rather than silently recomputing them with
  a possibly different root.
- For milestone-5 dual summaries, it is acceptable to use compact symbolic
  lines such as `phi_hat o phi = [n]` and short `yes` / `no` verdict lines, as
  long as the surrounding text makes clear that the checks were exhaustive on
  enumerated rational points.
- Avoid brittle full-output formatting tricks unless the exact layout is part
  of the teaching goal.

## Testing expectations

- Test important phrases and mathematical content, not giant exact snapshots,
  unless the output format is intentionally fixed.
- For explanatory helpers, test both the “works” case and the “does not exist /
  not supported / no exact answer” case when applicable.
- Keep tests aligned with the actual backend semantics, especially for
  approximate complex arithmetic and exact-only rational helpers.
- For isogeny helpers, test both structural summaries and at least one
  concrete small-field explanation of how a codomain or image point is
  computed.
- For milestone-7 helpers, test both the compact summary surface and the
  richer torsion explanation surface for the phrases the runnable example
  relies on.
- When a helper is intended to support a runnable example, test the important
  phrases that the example relies on rather than snapshotting the entire final
  console output.

## Review heuristics

A good change under `src/visualization` should improve at least one of:

- readability
- mathematical honesty
- pedagogical usefulness
- consistency with the underlying domain API

If a visualization helper hides an important caveat from the reader, it is not
done yet.
