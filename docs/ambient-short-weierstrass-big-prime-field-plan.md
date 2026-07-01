# Ambient Short-Weierstrass Over `BigPrimeField`

## Goal

Add one runtime-owned short-Weierstrass curve lane over
`fields::BigPrimeField` without refactoring the existing static
`F: Field` curve stack.

The milestone should make it possible to:

- build a nonsingular curve over one named large prime such as
  `2^255 - 19` or the secp256k1 prime,
- validate affine points,
- lift `x`-coordinates through the runtime square-root backend,
- add, negate, double, and scalar-multiply affine points,
- run at least one end-to-end educational example over a large prime.

Non-goal:

- do not try to make these runtime curves automatically implement the current
  `CurveModel`, `GroupCurveModel`, `FiniteGroupCurveModel`, or
  `FrobeniusTraceCurveModel` traits,
- do not attempt full parity with the static curve ecosystem in this
  milestone,
- do not introduce Schoof, Mestre, exhaustive enumeration, or isogeny
  plumbing for the runtime lane yet.

## Why A Separate Lane

The current curve stack is statically typed around `F: Field`:

- `CurveModel::BaseField` is constrained by `Field`,
- `AffinePoint<F>` stores `F::Elem`,
- concrete curve families such as `ShortWeierstrassCurve<F>` are parameterized
  by the field family type itself.

`BigPrimeField` is a runtime ambient object implementing `AmbientField`, not
`Field`, so it does not fit that stack automatically.

Trying to generalize the whole curve hierarchy immediately would create a much
larger and riskier refactor than the educational value of the first runtime
large-prime curve milestone justifies.

## Deliverable

Create a new runtime curve module, tentatively:

`src/elliptic_curves/ambient_short_weierstrass/`

with a small parallel surface:

- `AmbientShortWeierstrassCurve`
- `AmbientShortWeierstrassPoint`
- runtime group-law helpers
- focused tests
- one runnable example

This module should be explicit that it is a runtime-owned lane, not a drop-in
instance of the existing static trait family.

## Cargo Feature Gate

Compile this milestone behind the `runtime-field-curves` feature.

This feature owns the runtime large-prime field surface and the ambient
short-Weierstrass curve lane built on top of it. The name is intentionally a
little broader than `big-prime-field`: the field backend is useful on its own,
but its main purpose in this milestone is to make curves over runtime-selected
large primes possible.

Expected commands:

```text
cargo test -q --features runtime-field-curves big_prime_field
cargo run -q --features runtime-field-curves --example big_prime_field
```

Future examples and tests for `ambient_short_weierstrass` should use the same
feature gate.

## Proposed File Layout

```text
src/elliptic_curves/ambient_short_weierstrass/
  mod.rs
  curve.rs
  point.rs
  membership.rs
  lift_x.rs
  group_law.rs
  scalar_mul.rs
  display.rs
  tests.rs
examples/
  ambient_big_prime_short_weierstrass.rs
```

If the testing surface grows, later split `tests.rs` into:

- `tests/construction.rs`
- `tests/lift_x.rs`
- `tests/group_law.rs`
- `tests/large_prime_examples.rs`

## Data Model

### Curve

```rust
pub struct AmbientShortWeierstrassCurve {
    field: BigPrimeField,
    a: BigPrimeFieldElem,
    b: BigPrimeFieldElem,
}
```

Rules:

- `a` and `b` must already be elements of `field`,
- construction validates characteristic different from `2` and `3`,
- construction rejects singular curves via `Δ = -16(4a^3 + 27b^2) = 0`.

The field should be stored by value, not by reference, so curve values stay
self-contained and easy to pass around in examples and tests.

### Point

```rust
pub enum AmbientShortWeierstrassPoint {
    Infinity,
    Finite {
        x: BigPrimeFieldElem,
        y: BigPrimeFieldElem,
    },
}
```

Keep the explicit infinity variant, matching the educational style of the
static affine point layer.

Do not reuse `AffinePoint<F>` in this milestone: it is tied to `F: Field`.

## Milestone Phases

### Phase 1: Descriptor And Membership

Implement:

- `AmbientShortWeierstrassCurve::new(field, a, b)`,
- `field()`, `a()`, `b()`,
- discriminant,
- `contains(point)`,
- `identity()`,
- `is_identity(point)`,
- `to_equation_string()`.

Verification:

- singular rejection,
- characteristic `2` and `3` rejection,
- valid small-prime sample curve acceptance,
- point membership tests over a small runtime prime.

### Phase 2: Affine Point Construction And `x`-Lifting

Implement:

- `point(x, y) -> Result<AmbientShortWeierstrassPoint, CurveError-like error>`,
- `rhs_value(x) = x^3 + ax + b`,
- `lift_x(x)`,
- `point_from_x(x)`,
- `sqrt_pair`-driven lifting through `BigPrimeField`.

Notes:

- keep the same semantic shape as the static `LiftXCoordinate` story,
- but do not force the runtime lane into the existing trait.

Verification:

- one and two-point fibers over `F_17`,
- no-point cases,
- `y = 0` tangency case,
- one large-prime sanity case such as `x = 0` on a hand-picked curve.

### Phase 3: Native Affine Group Law

Implement the classical short-Weierstrass affine formulas:

- negation `(x, y) ↦ (x, -y)`,
- distinct-point addition,
- tangent doubling,
- identity handling,
- inverse-pair cancellation.

Keep the formulas local and fully documented in rustdocs near the executable
helpers.

Verification:

- hand-checked examples over `F_17`,
- identity and inverse laws,
- doubling of a `2`-torsion point returns infinity,
- associativity checked on a tiny fixed corpus over a small prime.

### Phase 4: Scalar Multiplication

Implement:

- `mul_scalar(point, u64)`,
- optionally `mul_scalar_biguint(point, &BigUint)` if it helps later large-prime
  workflows.

Recommendation:

- land `u64` first because it matches current educational curve APIs,
- add `BigUint` immediately only if the implementation stays small and clearly
  useful for the example surface.

Verification:

- repeated-addition agreement on tiny curves,
- `[0]P = O`,
- `[1]P = P`,
- `[n+m]P = [n]P + [m]P` on a small regression corpus.

### Phase 5: Educational Example

Add one runnable example:

`examples/ambient_big_prime_short_weierstrass.rs`

Suggested story:

1. construct `F_p` for one named large prime,
2. build one simple nonsingular curve,
3. try a few `x`-lifts,
4. choose one point,
5. compute `[n]P`,
6. print the curve equation, the chosen point, and the result.

Good defaults:

- a small demonstrative curve over the secp256k1 prime or `2^255 - 19`,
- a scalar that still prints quickly and deterministically.

## Error Surface

Do not overload the current static-field `CurveError` unless the fit is honest.

Two acceptable options:

1. add a small new runtime-curve error type under the ambient module,
2. reuse `CurveError` only if each runtime failure maps naturally to an
   existing variant.

My recommendation:

- start with a small local `AmbientShortWeierstrassError`,
- add `From<FieldError>` where appropriate,
- degrade into broader `CurveError` only later if a shared abstraction becomes
  genuinely valuable.

## Reuse Opportunities

Reuse from the current repo:

- mathematical formulas and testing ideas from static short-Weierstrass,
- `BigPrimeField` arithmetic, `sqrt`, and quadratic-character support,
- display tone and example style from existing `examples/`.

Avoid reusing too aggressively:

- `AffinePoint<F>`,
- `CurveModel` and `GroupCurveModel`,
- finite-field counting or enumeration traits.

## Verification Plan

Prefer focused verification at each phase:

- `cargo test -q ambient_short_weierstrass`
- if split by test names, use the narrowest honest filter available
- run the example directly once the example phase lands

Minimum final evidence for the milestone:

- construction tests pass,
- `lift_x` tests pass,
- group law tests pass,
- scalar multiplication tests pass,
- example runs successfully.

## Risks

### Risk 1: Over-generalizing too early

Trying to reuse the static trait stack too early may cause more churn than the
milestone value justifies.

Mitigation:

- keep the runtime lane separate in this milestone.

### Risk 2: Duplicated formulas diverge from static short-Weierstrass

Once two executable short-Weierstrass lanes exist, formulas can drift.

Mitigation:

- document formulas explicitly,
- keep tests over small primes exhaustive or near-exhaustive where feasible,
- consider later extracting shared formula helpers only after the runtime lane
  stabilizes.

### Risk 3: Error-surface confusion

Mixing `FieldError`, `CurveError`, and runtime-lane-specific failures without a
clear boundary will make the API harder to read.

Mitigation:

- pick one local runtime error type early,
- wire conversions intentionally.
  g

## Out Of Scope For This Milestone

- runtime Montgomery or general-Weierstrass families,
- reductions between runtime curve families,
- projective coordinates,
- Hasse interval search,
- point order and group order algorithms,
- isogenies,
- integration with the static generic curve traits.

## Follow-up Milestones

If this milestone succeeds, likely next steps are:

1. `Ambient short-Weierstrass + BigUint scalar multiplication`
2. `Runtime point-order helpers from known multiples`
3. `Runtime Montgomery reduction or same-field conversion stories`
4. evaluate whether a small `AmbientCurveModel` trait family is worth adding

## Recommendation

Proceed with this milestone only as a separate runtime lane.

That gives the repo real value:

- curves over named large primes,
- end-to-end arithmetic,
- educational examples that exercise `BigPrimeField`,

without paying the much larger cost of re-architecting the entire existing
curve stack in the same pass.
