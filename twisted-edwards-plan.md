# Twisted Edwards Introduction Roadmap

This document records the current staged plan for adding
`TwistedEdwardsCurve<F>` to the repo.

The mathematical target is the affine family

`a x^2 + y^2 = 1 + d x^2 y^2`

initially over fields of characteristic different from `2`.

The goal is not to land the whole Edwards ecosystem at once. The goal is to
introduce the model in slices that preserve the repo's current strengths:

- explicit model descriptors
- honest capability boundaries
- reusable shared traits where they fit
- explicit cross-model conversion witnesses
- small-field verification before deeper algorithm layers

## Repository state today

The existing curve-model stack already gives us several strong anchors:

- `ShortWeierstrassCurve<F>` is the mature executable family.
- `GeneralWeierstrassCurve<F>` already follows the staged-model pattern:
  validated descriptor, invariants, membership, `LiftXCoordinate`, native
  group law, projective layer, and explicit reductions.
- `MontgomeryCurve<F>` already owns:
  - validated descriptors for `B y^2 = x^3 + A x^2 + x`
  - invariants and membership
  - `LiftXCoordinate`
  - native affine group law
  - explicit conversion to short and general Weierstrass
  - small finite-field wrappers and an educational example
- cross-model conversion infrastructure currently assumes a reusable witness
  with total point transport in both directions:
  `CurveModelConversion`

That makes Montgomery the natural bridge for Twisted Edwards at the
coefficient level, because the classical formulas are direct:

- from Twisted Edwards to Montgomery:
  - `A = 2(a + d) / (a - d)`
  - `B = 4 / (a - d)`
- from Montgomery to Twisted Edwards:
  - `a = (A + 2) / B`
  - `d = (A - 2) / B`

So whole-curve conversion is algebraically cheap and does not require square
roots or field enumeration.

## Adopted design decisions

This pass closes the main architectural questions for the first milestone.

### Decision 1: identity is the finite affine point `(0, 1)`

Adopted policy:

- keep `TwistedEdwardsCurve<F>` on the shared `AffinePoint<F>` value type
- represent the neutral element honestly as the finite point `(0, 1)`
- do not invent an artificial point at infinity just to mimic Weierstrass

Why:

- it matches the natural affine geometry of the model
- it keeps the repo educationally honest
- it avoids the anti-pattern of forcing every curve family into a hidden
  Weierstrass semantic mold

Immediate consequence:

- shared finite enumeration must become identity-aware before Twisted Edwards
  leans on blanket `EnumerableCurveModel` behavior
- for Twisted Edwards, the conceptual point set is "all affine solutions",
  with the model itself knowing which one is the identity

### Decision 2: curve conversion is total, point transport starts partial

Adopted policy:

- whole-curve Twisted-Edwards ↔ Montgomery conversion is a first-milestone
  feature
- point transport starts as honest birational transport on the affine open
  subset where the classical formulas are defined
- total rational-point correspondence is a later milestone that will need a
  richer abstraction or a different point representation

Why:

- coefficient conversion is clean and global
- point conversion is birational, with real exceptional loci
- pretending a global affine point equivalence exists would make the API lie

### Decision 3: first native group law stays generic in `(a, d)`

Adopted policy:

- the first executable group law is for generic valid twisted-Edwards
  coefficients `(a, d)`
- the repo does not specialize early to `a = -1`
- the first affine formulas are documented as honest formulas with possible
  denominator failure, not as "complete" formulas

Why:

- it preserves the educational value of the full family
- it avoids overfitting the first milestone to one cryptographic subfamily
- it keeps later "complete formulas" or special-subfamily work clearly
  separated

## Three levels that must stay distinct

This distinction is now part of the intended design contract:

1. Curve conversion
   Coefficient-level transport between model descriptors. For Twisted Edwards
   and Montgomery, this should be total in the first milestone.
2. Birational point transport on an affine open
   Point transport using the classical rational formulas only where they are
   actually defined. This should be explicit and partial.
3. Total rational-point correspondence
   A future, richer story that accounts for exceptional points honestly and
   should not be conflated with the previous level.

Keeping these layers separate now should prevent later confusion around
Montgomery ladders, torsion, cofactors, Decaf/Ristretto-style quotients, and
isogeny-side transport.

## What should work "for free" and what should not

Once the baseline model exists, the following should be easy or nearly free:

- `CurveModel`
- `AffineCurveModel`
- `HasJInvariant`
- whole-curve Montgomery conversion in both directions
- whole-curve conversion to short/general Weierstrass by composition
- `LiftXCoordinate` from the fiber equation
  `y^2 = (1 - a x^2) / (1 - d x^2)`
- exhaustive small-field enumeration once the shared enumerable path handles a
  finite identity honestly

The following are not free and should be staged explicitly:

- total point transport Edwards ↔ Montgomery
- complete Edwards addition formulas without coefficient restrictions
- projective or extended Edwards coordinates
- Edwards-native finite-field order algorithms beyond exhaustive routes
- Edwards-native isogeny, function-field, or division-polynomial layers

## Recommended sequence

## Stage 0: shared-trait fit check

### Why first

Twisted Edwards is the first obvious affine family in this repo whose neutral
element is naturally finite. That makes it a useful pressure test for the
shared affine/enumeration traits.

### Deliverables

- confirm that `CurveModel`, `AffineCurveModel`, `LiftXCoordinate`, and
  `GroupCurveModel` can support a finite identity without semantic drift
- adjust shared enumerable behavior so identity points lifted from an affine
  fiber are not duplicated in `points()`
- document the chosen point-representation story for the new family

### Exit criteria

- the repo has one explicit documented policy for affine models with finite
  identity
- no shared helper silently assumes "identity = infinity"

## Stage A: validated curve descriptor and invariants

### Deliverables

- add `models/twisted_edwards/type_definition.rs`
- define `TwistedEdwardsCurve<F>` for
  `a x^2 + y^2 = 1 + d x^2 y^2`
- validate the first supported regime:
  - characteristic different from `2`
  - `a != 0`
  - `d != 0`
  - `a != d`
- expose coefficient getters and `to_equation_string()`
- add classical invariants or compute them through the direct Montgomery
  companion formulas when that keeps the code smaller and still explicit

### Notes

For this family, the non-singularity conditions are cheap and explicit, so
this stage should stay very small.

### Exit criteria

- construction rejects singular/unsupported inputs honestly
- one small test file covers construction and invariant sanity checks

## Stage B: membership, model traits, and `lift_x`

### Deliverables

- add equation-membership helpers
- implement:
  - `CurveModel`
  - `AffineCurveModel`
  - `HasJInvariant`
  - `LiftXCoordinate`
- use the direct fiber equation
  `y^2 = (1 - a x^2) / (1 - d x^2)`
- treat the vanishing denominator honestly rather than silently cancelling it

### Why this is a good early win

This unlocks enumeration and visualization immediately, and it reuses the same
pedagogical `x`-fiber story already used by Montgomery and general
Weierstrass.

### Exit criteria

- `point`, `contains`, and `point_from_x` work over supported `SqrtField`
  backends
- tiny-field tests cover both generic fibers and the neutral point `(0, 1)`

## Stage C: whole-curve Montgomery conversion

### Deliverables

- add a model-owned conversion layer, likely `reduction.rs`
- expose:
  - `as_montgomery()` from Twisted Edwards
  - `as_twisted_edwards()` from Montgomery
  - `From`/`TryFrom` impls for whole curves in both directions
- document the coefficient formulas explicitly
- keep the witness reusable for later composed whole-curve conversions to
  short and general Weierstrass

### Important scope boundary

This stage is about curve descriptors first, not about total point transport.

### Exit criteria

- whole-curve roundtrips preserve coefficients
- composed conversions to short/general Weierstrass agree with direct
  invariant checks

## Stage D: birational point transport on the affine open

### Deliverables

- add an explicit partial point-transport API
- make the partial nature clear in naming and docs
- document the standard formulas and their excluded loci
- add roundtrip tests only on the certified domain of definition

### Recommended API direction

Prefer names that make the semantics visible, for example:

- `try_point_to_montgomery_open(...)`
- `try_point_from_montgomery_open(...)`

or a similarly explicit birational naming surface.

### Recommended first formulas

For finite Edwards points away from `y = 1` and `x = 0`, use the classical
bridge

- `u = (1 + y) / (1 - y)`
- `v = (1 + y) / (x (1 - y))`

and for Montgomery points away from `u = -1` and `v = 0`, use

- `x = u / v`
- `y = (u - 1) / (u + 1)`

but keep the API honest about exceptional points instead of pretending these
formulas define a total affine equivalence.

### Exit criteria

- the repo has one explicit documented Montgomery/Edwards birational
  point-transport policy
- tests show both successful transport and honest rejection on exceptional
  inputs

## Stage E: native affine group law

### Deliverables

- add a native affine `GroupCurveModel` implementation for Twisted Edwards
- document the formulas used and their preconditions
- keep denominator-failure behavior explicit in the first implementation

### Recommended first formulas

Use the generic affine formulas

- `x3 = (x1*y2 + y1*x2) / (1 + d*x1*x2*y1*y2)`
- `y3 = (y1*y2 - a*x1*x2) / (1 - d*x1*x2*y1*y2)`

and do not describe them as complete unless a later restricted subfamily is
introduced with that guarantee documented explicitly.

### Validation strategy

Prefer native formulas for the executable path. Use Montgomery transport only
as a test oracle where the birational domain allows it, not as the permanent
execution core.

### Exit criteria

- `neg`, `add`, `double`, and `mul_scalar` work on supported coefficients
- tests cover identity, inverse, doubling, and small-field group laws

## Stage F: finite-field integration and shared wrappers

### Deliverables

- ensure `EnumerableCurveModel` now works honestly with finite identity
- let small-field routes come online through shared traits:
  - `points()`
  - `order()`
  - `point_order()`
  - `group_structure()`
  - `group_exponent()`
- decide which higher-level wrappers should be model-owned immediately and
  which ones can wait

### Recommendation

Match the current general/Montgomery posture:

- exhaustive small-field routes can be native immediately
- deeper routes should stay deferred until the model owns enough algebra to
  justify them

### Exit criteria

- dedicated compatibility tests confirm that shared finite-group helpers work
  with a finite affine identity

## Stage G: examples, visualization, and public-root wiring

### Deliverables

- add one educational example analogous to `examples/montgomery.rs`
- add visualization helpers for:
  - the Twisted Edwards equation
  - the neutral point and inverse symmetry
  - the Montgomery companion
  - one transported sample calculation
- wire the new model into:
  - `src/elliptic_curves/models/mod.rs`
  - `src/elliptic_curves/mod.rs`

### Exit criteria

- the model is discoverable from the crate root
- one runnable example demonstrates both native Edwards behavior and the
  Montgomery bridge

## Later milestones that should stay separate

- total rational-point correspondence across exceptional loci
- projective or extended Edwards coordinates
- complete addition formulas under documented coefficient restrictions
- x-only ladders or Edwards-side scalar schedules
- Edwards-native reduction stories beyond Montgomery composition
- Edwards-native function fields, isogenies, or division polynomials
