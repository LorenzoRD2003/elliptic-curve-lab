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
  - `GroupCurveModel` for models that expose actual point addition, doubling,
    and scalar multiplication
  - `LiftXCoordinate` for models that can recover points from `x`
  - `EnumerableCurveModel` only for small finite settings where exhaustive
    point listing is honest
  - `FiniteGroupCurveModel` only for small finite settings where point orders
    and related order-theoretic helpers can be computed by direct group traversal
  - invariant capability traits such as `HasJInvariant` when a family can
    expose classical invariants honestly without inflating `CurveModel`
- A small `CurveIsomorphism` trait is now part of the intended architecture for
  explicit curve-to-curve witnesses. It should stay narrow: domain, codomain,
  and point evaluation.
- Exhaustive base-field isomorphism search on enumerable finite fields is now
  a first-class educational tool, not just an internal convenience. It is
  acceptable to use it to support milestone-5 workflows such as dual-isogeny
  search, provided the docs say clearly that this is a tiny-field exhaustive
  routine.

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
- When an invariant is useful across multiple consumers but still not universal,
  prefer a narrow capability trait over adding it directly to `CurveModel`.
- Short-Weierstrass coefficient-scaling helpers such as `scaled_by` and
  `isomorphic_via_scale` belong on `ShortWeierstrassCurve` itself, since they
  describe model-specific coefficient transport rather than a generic
  curve-morphism interface.
- Explicit short-Weierstrass scaling witnesses belong in
  `ShortWeierstrassIsomorphism<F>` and currently cache their codomain curve as
  validated state. That means the `CurveIsomorphism` trait may return the
  codomain by reference honestly, without recomputing it.
- For short-Weierstrass isomorphisms, use the convention
  `\phi_u : E -> E'`, `(x, y) -> (u^2 x, u^3 y)`.
  If `E : y^2 = x^3 + ax + b`, then document and implement the image model as
  `E' : y^2 = x^3 + a'x + b'` with `a' = u^4 a` and `b' = u^6 b`.
  Treat this as the canonical normalization for M4 unless a later milestone
  explicitly introduces a second convention and explains the translation.
- Keep the distinction explicit between:
  same `j`-invariant = isomorphic over an algebraic closure,
  versus
  isomorphic over the base field = there exists `u in F^*` with
  `a' = u^4 a` and `b' = u^6 b`.
  Do not collapse those two notions in docs or API names.
- For quadratic-twist objects, store the twist factor `d` as the primary data.
  Do not store a base-field scaling witness `u` as mandatory state, since that
  witness may fail to exist exactly in the genuinely quadratic case.
- Point enumeration is acceptable only when the base field is explicitly small
  and enumerable. Say so in docs.
- For `EnumerableFiniteField`, exhaustive witness search is acceptable for
  pedagogical helpers such as finding a concrete short-Weierstrass scaling
  isomorphism or enumerating all compatible base-field isomorphisms between two
  short-Weierstrass models. Keep the docs honest that this is a small-field
  educational routine, not a large-field optimized algorithm.
- Point-order and torsion helpers are acceptable only when the ambient group is
  explicitly small and enumerable. Say directly that the current algorithms
  use direct traversal or repeated addition rather than efficient large-group
  techniques.
- Do not rush into optimized formulas, scalar multiplication, serialization, or
  cryptographic hardening.
- If a model exposes a group-law trait, keep the docs explicit about which
  operations are baseline educational formulas and which errors represent
  invalid off-curve inputs.
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
- When group operations are exposed, test identity, inverses, doubling, scalar
  multiplication, and at least one small exact associativity example.
- When point-order or torsion helpers are exposed, test at least one identity
  case, one non-trivial finite-order example, and one invalid off-curve input.
- When a helper depends on field-side capabilities, add at least one test that
  exercises the positive path and one that shows the honest negative path.
- For enumeration helpers, test the identity case, finite-point count, and at
  least one small exact order example.
- For short-Weierstrass isomorphism comparisons, include explicit tests for
  the special `j = 0` (`a = 0`) and `j = 1728` (`b = 0`) families, in
  addition to generic `a,b != 0` examples.
- When short-Weierstrass isomorphisms are cached objects rather than
  recomputed witnesses, test both the coefficient transport and the cached
  codomain-facing behavior through the generic `CurveIsomorphism` trait
  surface.
- For short-Weierstrass automorphism helpers over enumerable finite fields,
  test the generic case separately from the `j = 1728` and `j = 0` special
  families, since those special loci can admit extra automorphisms.
- For quadratic-twist helpers, test at least:
  preservation of the `j`-invariant,
  one square-factor case that stays base-field isomorphic,
  one non-square case that is not base-field isomorphic in the chosen sample,
  and the point-count relation `#E(F_p) + #E^(d)(F_p) = 2p + 2` over small odd
  prime fields.

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
