# Montgomery Ladder Roadmap

This document starts after the first `MontgomeryCurve<F>` milestone is already
complete.

That baseline now exists:

- validated descriptors for `B y^2 = x^3 + A x^2 + x` in characteristic
  different from `2`
- classical invariants
- direct `lift_x`
- explicit conversion to short and general Weierstrass companions
- native affine group law
- small finite-field wrappers, compatibility tests, visualization, and one
  runnable example

The next stage is no longer “make Montgomery exist”. It is “make Montgomery
own its distinctive `x`-coordinate arithmetic story”.

The main educational target is:

- compute `x([n]P)` from `x(P)` without requiring the sign of `y(P)`
- make the differential-addition structure explicit
- prepare the repo for later stories around Curve25519-style APIs,
  side-channel-aware scalar-multiplication discussions, and efficient
  Montgomery-model workflows

This roadmap is intentionally dependency-ordered. The goal is to reach a real
`MontgomeryLadder` surface without smuggling in half-verified formulas or
collapsing everything back into affine addition.

## Guiding rules

1. Keep the ladder story separate from the affine group law.
   The point of this milestone is not “yet another scalar multiplication
   wrapper”; it is the `x`-only differential arithmetic itself.
2. Preserve mathematical honesty about inputs and outputs.
   If a routine consumes only `x(P)`, say so. If it additionally needs
   `x(P-Q)` or a checked affine witness, say so too.
3. Keep the repo’s general `B y^2 = x^3 + A x^2 + x` model explicit.
   Classical ladder formulas are often presented for `B = 1`. If an internal
   normalization or rescaling is used, document exactly what is and is not
   available over a general base field.
4. Treat “ladder-shaped” and “constant-time” as different claims.
   The first milestone should aim for a fixed-operation educational ladder. A
   stronger constant-time claim belongs only after the implementation details
   and API surface honestly support it.
5. Close each stage with cross-model tests.
   Whenever an `x`-only formula is added, compare it against affine
   scalar multiplication whenever a full affine point is available.

## What is still intentionally deferred

The following items should not be bundled into the first ladder milestone:

- full Edwards/Montgomery cryptographic API design
- byte encoding / decoding formats
- clamped scalar APIs
- production-grade constant-time guarantees
- Curve25519-specific field backends or protocol code
- Montgomery-only isogeny or function-field work

## Recommended sequence

## Stage A: `x`-coordinate value objects and ownership boundaries

### Why first

Right now the Montgomery model owns affine points and affine arithmetic, but it
does not yet own a first-class `x`-only representation. The ladder should not
be introduced as a loose helper that passes raw field elements around without
semantics.

### Deliverables

- introduce an explicit `x`-coordinate representation for Montgomery work,
  likely one affine-facing type for a checked finite `x` value and one
  projective `X:Z` representation for ladder execution
- decide how the point at infinity is represented on the `x`-line and in the
  projective ladder state
- document whether the value object means:
  - an affine `x` of some finite point,
  - a projective `X:Z` ratio up to scale, or
  - an `x`-line class that forgets the sign of `y`
- keep this layer model-owned under `models/montgomery/` rather than as a new
  generic top-level abstraction too early
- pair that `x`-coordinate layer with the explicit normalization witness above,
  so the eventual ladder execution can live on
  `NormalizedMontgomeryCurve<F>` while the outer `MontgomeryCurve<F>` keeps the
  honest “available only when `B` is a square” wrapper

### Design fork to resolve here

The repo’s public Montgomery model is `B y^2 = x^3 + A x^2 + x`, while many
classical ladder formulas are written for `v^2 = u^3 + A u^2 + u`.

At this stage, decide which execution story the repo will own first:

- native `B`-aware `X:Z` formulas, or
- an explicit internal normalization story, only when the required scaling is
  mathematically available over the current base field

Do not blur those two stories together.

### Chosen direction

The first ladder milestone should use the second route:

- an explicit internal normalization story
- only when the required scaling really exists over the current base field
- with honest failure or unavailability when that normalization witness cannot
  be produced

This keeps the first `x`-only implementation close to the classical Montgomery
ladder formulas while preserving mathematical honesty about which base fields
actually support that normalization.

More concretely, the intended normalization witness is the `y`-rescaling
`v = sqrt(B) y`, which rewrites
`B y^2 = x^3 + A x^2 + x`
as the normalized Montgomery model
`v^2 = x^3 + A x^2 + x`
over the same base field.

Therefore:

- the normalization exists over `K` if and only if `B` is a square in `K`
- if `B` is not a square in `K`, the source model should be treated honestly as
  a quadratic twist of the normalized Montgomery model over `K`
- the first ladder API should not hide that twist by silently pretending the
  normalized formulas apply without a normalization witness

One good target shape for this stage is:

```rust
pub struct NormalizedMontgomeryCurve<F> {
    a: F::Elem,
}

pub struct MontgomeryNormalization<F: Field> {
    source: MontgomeryCurve<F>,
    target: NormalizedMontgomeryCurve<F>,
    sqrt_b: F::Elem,
}

impl<F: Field + SqrtField> MontgomeryCurve<F> {
    pub fn try_normalize(
        &self,
    ) -> Result<MontgomeryNormalization<F>, MontgomeryNormalizationError>;
}
```

The final error surface does not need to be fixed yet. It may be one dedicated
Montgomery-local error or one Montgomery-specific branch integrated honestly
with the repo’s broader `CurveError` surface.

### Exit criteria

- one explicit `x`-coordinate type exists
- one explicit normalization witness exists for the `B`-square case
- the docs say what information that type remembers and what it forgets
- tests cover basic construction, equality-up-to-scale, and infinity semantics
- tests also cover successful normalization when `B` is a square and honest
  rejection when `B` is not a square over the current base field

## Stage B: Differential arithmetic primitives

### Why here

The ladder is powered by differential formulas, not by generic affine
addition. Before introducing any bit-walking algorithm, the repo should own the
three core ingredients as named, tested operations.

### Deliverables

- add native Montgomery `xDBL`
- add native Montgomery differential addition `xADD`, with an explicit precondition
  that the difference point is known through its `x`-coordinate data
- optionally add a combined `xDBLADD` helper if that is the clearest execution
  unit for the eventual ladder
- document the exact formulas used, the coordinate chart, and the required
  preconditions
- if one cached coefficient such as an `A24`-style term is introduced, make it
  a named Montgomery-owned helper with docs explaining why it appears

### Test strategy

- compare `xDBL` against the affine `double` result when a full affine point is
  known
- compare `xADD` against affine `P + Q` on examples where `P`, `Q`, and
  `P - Q` are all known
- include tiny-field exhaustive tests over representative supported
  characteristics such as `3` and `5`

### Exit criteria

- `xDBL` and `xADD` are available as documented model-owned operations
- the implementation is validated against affine arithmetic wherever affine
  witnesses exist
- no ladder code exists yet that duplicates these formulas ad hoc

## Stage C: Ladder state machine and scalar schedule

### Why here

Once `xDBL` and `xADD` are stable, the actual ladder becomes a small state
machine rather than a formula dump.

### Deliverables

- define a `MontgomeryLadder` execution surface or equivalent model-owned
  helper layer
- represent the standard invariant pair, e.g. neighboring multiples whose
  difference is the original base point
- implement the bit-by-bit ladder schedule using only the differential
  primitives
- keep the first scalar type simple and explicit, for example `u64` or one
  existing repo integer surface, before generalizing to larger scalars
- expose the core educational output:
  `x([n]P)` from `x(P)`

### API candidates

One good first milestone would be an API family like:

- `ladder_x(&self, base_x: ..., scalar: u64) -> Result<..., CurveError>`
- `ladder_x_with_trace(...) -> LadderReport<...>`
- an internal helper that returns the final ladder pair, not just the final
  `x([n]P)`, if that makes verification and teaching easier

If Stage A keeps the normalization witness explicit, a natural wrapper shape is:

```rust
impl<F: Field + SqrtField> MontgomeryCurve<F> {
    pub fn try_ladder_x(
        &self,
        x: F::Elem,
        n: u64,
    ) -> Result<F::Elem, MontgomeryNormalizationError> {
        let normalization = self.try_normalize()?;
        Ok(normalization.target.ladder_x(x, n))
    }
}
```

because the `x`-coordinate is unchanged by the `y`-rescaling normalization.

That does not force the final public API to use this exact spelling, but it is
the intended architectural direction for the first milestone.

### Important honesty rule

If the implementation is operation-regular but not yet audited as constant-time
with respect to the host language and field backend, document it as
ladder-shaped or fixed-schedule, not as production constant-time.

### Exit criteria

- the repo can compute `x([n]P)` from `x(P)` without a `y` sign
- the ladder reuses only the Stage B primitives
- tests compare ladder output against affine scalar multiplication whenever a
  full point lift is available

## Stage D: Validation, reports, and characteristic-sensitive coverage

### Why here

The ladder is subtle because it forgets information. The educational story gets
much stronger if the repo shows exactly what can still be certified and what
cannot.

### Deliverables

- add a focused `tests/ladder.rs` or similar suite under Montgomery
- add examples where two affine points `P` and `-P` share the same input
  `x(P)` and therefore the same ladder output `x([n]P)`
- add property tests comparing:
  - `ladder_x(x(P), n)` and `x([n]P)` from affine scalar multiplication
  - `ladder_x(x(P), n)` and `ladder_x(x(-P), n)`
- add docs or a small report type explaining when `x([n]P)` is enough and when
  recovering the full point still needs extra data

### Exit criteria

- the limitation “the ladder returns an `x`-coordinate class, not a signed
  affine point” is explicit in docs and tests
- representative tiny-field exhaustive checks and broader property coverage are
  in place

## Stage E: Ergonomic public wrappers and visualization

### Why here

Once the executable ladder is trustworthy, the next step is making it visible
and educational for callers.

### Deliverables

- add ergonomic curve-side wrappers on `MontgomeryCurve<F>`
- add visualization helpers that show:
  - the source Montgomery equation
  - the input `x(P)`
  - the scalar
  - the resulting `x([n]P)`
  - optionally the final ladder pair or the differential-invariant story
- add one runnable example centered on `x([n]P)` rather than generic affine
  addition
- keep docs explicit about which parts are native Montgomery `x`-arithmetic and
  which comparisons are only for validation through affine points

### Exit criteria

- one example demonstrates the `x`-only scalar story end to end
- visualization/docs explain the differential-addition invariant in plain
  language

## Stage F: Side-channel-aware educational layer

### Why here

This repo is educational, so there is real value in making the “why the ladder
matters” story first-class, even before any production cryptography surface is
attempted.

### Deliverables

- add a short educational note explaining why ladder schedules are attractive
  for regular scalar multiplication
- distinguish carefully between:
  - regular control flow
  - constant-time aspirations
  - backend-dependent leakage that the repo does not yet control
- if appropriate, add a traced execution report that illustrates the repeated
  `xDBL` / `xADD` pattern without claiming security properties the code has not
  earned

### Exit criteria

- the docs connect the ladder to side-channel-aware scalar-multiplication
  design without overstating guarantees

## Stage G: Future bridges

### Why last

Only after the core ladder exists should the repo decide how far to lean into
Montgomery-specific cryptographic workflows.

### Possible next steps after the ladder milestone

- a canonical “Montgomery `u`-coordinate” educational wrapper for models where
  that naming is clearer
- larger-scalar support beyond `u64`
- optional y-recovery workflows when extra data is available
- transport/comparison examples against short-Weierstrass scalar multiplication
- a future Curve25519-oriented roadmap:
  field-specific normalization, scalar conventions, encodings, and protocol
  stories

### Exit criteria

- the core ladder milestone is complete before any protocol-shaped expansion is
  attempted

## Immediate recommendation

If work resumes immediately after this roadmap, the best next target is:

1. Stage A with an explicit Montgomery-owned `X:Z` value object.
2. Stage B with documented `xDBL` and `xADD`.
3. Stage C with one small, honest `ladder_x(..., scalar: u64)` surface.

That sequence gets you to the mathematically interesting part quickly while
still keeping the API and documentation honest.
