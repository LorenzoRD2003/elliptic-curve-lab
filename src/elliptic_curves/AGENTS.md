# AGENTS.md for `src/elliptic_curves`

## Module mission

The `elliptic_curves` module should introduce curve models gradually and
honestly.

Right now the goal is not to ship a production EC library. The goal is to make
curve structure, point representations, and equation checks easy to read and
easy to extend.

## Current posture

- Early-stage scaffolding is acceptable when it is explicit and tested.
- Short Weierstrass support is currently the main concrete path.
- The current affine representation should preserve mathematical invariants in
  the type when possible.
- Validation logic such as discriminant checks and point-membership checks is
  part of the educational API surface, not incidental glue.
- Classical short-Weierstrass invariants such as `Δ`, `c4`, `c6`, and
  `j` are appropriate here when their docs explain the mathematics directly.
- Curve-side capability traits are now part of the intended architecture:
  - `AffineCurveModel` for checked affine construction
  - `LiftXCoordinate` for models that can recover points from `x`
  - `EnumerableCurveModel` only for small finite settings where exhaustive
    point listing is honest

## Design priorities

- Mathematical honesty before feature count.
- Clear point representations before group-law optimization.
- Conservative public APIs that explain their preconditions.
- Small, verifiable steps.

## Representation rules

- Prefer representations that make invalid states hard to express.
- The point at infinity should be modeled explicitly rather than smuggled
  through meaningless affine coordinates.
- If a constructor claims to return a point on a curve, it should validate that
  claim.
- If a curve model is only valid away from special characteristics, document
  that fact directly in the type or constructor docs.

## Scope guidance

- It is fine to start with affine membership checks, discriminants, and simple
  point constructors.
- Model-specific invariants can stay as inherent methods when they belong only
  to one presentation, such as short-Weierstrass invariants on
  `ShortWeierstrassCurve`.
- Point enumeration is acceptable only when the base field is explicitly small
  and enumerable. Say so in docs.
- Do not rush into optimized formulas, scalar multiplication, serialization, or
  cryptographic hardening.
- If a new curve API depends on extra field capability, such as square roots,
  prefer a narrow trait bound like `SqrtField` over broadening unrelated base
  traits.

## Error conventions

- Keep recoverable curve-domain failures in `CurveError`.
- Prefer specific variants such as unsupported characteristic, singular curve,
  or point-not-on-curve over ad hoc strings.
- Add a new error variant only when it expresses a genuinely distinct curve
  failure mode.

## Testing expectations

- Test both valid and invalid curve construction.
- Test both valid and invalid point construction.
- Test the point at infinity behavior explicitly when it participates in the
  public model.
- When a helper depends on field-side capabilities, add at least one test that
  exercises the positive path and one that shows the honest negative path.
- For enumeration helpers, test the identity case, finite-point count, and at
  least one small exact order example.

## Documentation expectations

- Public curve items should explain the mathematical model they represent.
- If a formula is valid only in characteristic different from `2` and `3`, say
  so directly.
- If a feature is educational, partial, or not yet a full group-law layer, say
  so explicitly.
- If an invariant is attached to a specific curve presentation, document both
  its defining formula and its mathematical role.
- If a helper only makes sense for small finite fields, say so directly in the
  rustdocs.
- Use concrete examples where they clarify the model.

## Review heuristics

A good change under `src/elliptic_curves` should improve at least one of:

- invariant safety
- readability
- mathematical honesty
- test coverage

If a curve change makes the point model or equation semantics harder to
explain, it is probably moving too fast for the current phase.
