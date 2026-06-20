# General Weierstrass Deep-Generalization Roadmap

This document starts after the first `GeneralWeierstrassCurve<F>` milestone is
already complete.

That baseline now exists:

- validated curve descriptors
- classical invariants
- explicit conversion to and from a short-Weierstrass companion
- native affine group law
- `lift_x` through `y`-fiber solving
- small finite-field integration
- compatibility tests and one educational example

The next stage is no longer “make the model exist”. It is “decide which parts
of the stack should become honestly general-model-native, and in which order”.

This roadmap is intentionally dependency-ordered. The goal is to avoid
generalizing everything at once.

## Guiding rules

1. Preserve mathematical honesty.
   If a feature still works only through the short companion, say so in the
   API, docs, and tests.
2. Prefer native algebra only when it unlocks a real next layer.
   Do not rewrite a short-companion route just for symmetry.
3. Keep characteristic-specific behavior explicit.
   Characteristic `2`, characteristic `3`, and characteristic `> 3` should not
   silently collapse into one undocumented implementation story.
4. Ship one reusable layer at a time.
   Each stage should close with docs, tests, and at least one example or
   visualization update when the surface is user-facing.

## What is already intentionally deferred

The following items are still not first-class native general-model features:

- projective-coordinate group law
- general-model function fields
- general-model isomorphisms as an explicit surface
- native general-model isogenies
- division-polynomial tooling for the general model
- Schoof or other order algorithms written natively for the general model
- graph-layer witnesses owned by the general model

## Recommended sequence

## Stage A: Projective group law

### Why first

The current affine formulas are correct and well tested, but they are not the
long-term execution surface. Any serious future work on native isogenies,
division polynomials, or larger finite-field workflows will want a projective
story.

### Deliverables

- add a projective point representation for the general model, or a shared
  projective layer with model-specific formulas
- implement native projective addition, doubling, negation, and scalar
  multiplication for `GeneralWeierstrassCurve<F>`
- keep the current affine formulas as a checked reference path during the
  transition
- preserve the current explicit TODO explaining why projective formulas replace
  the affine ones

### Dependencies

- none beyond the current baseline

### Exit criteria

- public group operations use the projective route internally
- affine and projective results agree on exhaustive tiny-field tests
- characteristic `2`, `3`, and `> 3` remain covered

## Stage B: Explicit admissible isomorphisms of the general model

### Why here

Once the group law is no longer “temporary affine”, the next structural layer
should be the actual change-of-variables theory of general Weierstrass models.
That is the right abstraction boundary for later function fields, reductions,
and native isogeny normalization.

### Deliverables

- introduce a general-model isomorphism surface for admissible changes of
  variables
- make the current short-companion reduction visibly reuse that language
  instead of feeling like one isolated special case
- expose composition, inverse, and point transport
- document the exact coefficient transformation rules

### Dependencies

- Stage A is recommended but not strictly required

### Exit criteria

- short embedding and short reduction can be described as specific isomorphism
  witnesses or closely related admissible maps
- tests cover coefficient transport, point transport, identity, inverse, and
  composition

## Stage C: General-model function fields

### Why here

Function fields should sit on top of the model’s own equation and isomorphism
language, not directly on top of one special reduction story.

### Deliverables

- model `F(E)` for
  `y^2 + a1xy + a3y = x^3 + a2x^2 + a4x + a6`
- define multiplication using the general relation in a way that stays readable
- add substitution helpers tied to the general model
- keep the representation educational and explicit, as in the short model

### Dependencies

- Stage B strongly recommended

### Exit criteria

- values reduce correctly modulo the general equation
- the general function field admits point validation and basic arithmetic tests
- at least one example or visualization helper explains the basis and the
  relation used

## Stage D: Native general-model isogenies

### Why here

Native isogenies are meaningful only once the model has its own algebraic and
function-field language.

### Deliverables

- decide whether the first native surface should be:
  - explicit maps defined through admissible isomorphisms and transported short
    data, or
  - genuinely native formulas
- add a minimal, honest first layer instead of aiming immediately for full
  Vélu parity with short Weierstrass
- if the first step still transports through short companions, record that as a
  staged bridge, not as “native”

### Dependencies

- Stage B required
- Stage C strongly recommended for pullback and differential stories

### Exit criteria

- one nontrivial isogeny workflow can be expressed without ad hoc coordinate
  glue living outside the model
- tests distinguish transported routes from native ones

## Stage E: Division polynomials and torsion search

### Why here

This is the first major item where “generalize or keep short-specific?” needs a
real decision. Division-polynomial infrastructure is expensive to generalize,
so it should happen only if the general model already owns enough algebra to
justify it.

### Deliverables

- decide whether to define division-polynomial tooling:
  - directly on the general model, or
  - through a canonical reduction layer plus an explicit witness surface
- if native, define the exact polynomial/rational objects that belong to the
  general model
- reconnect exact-order torsion search, comparisons, and visual explanations

### Dependencies

- Stage B required
- Stage C recommended

### Exit criteria

- one honest torsion-search API exists for the general model
- documentation says whether the route is native or transported
- property tests compare against exhaustive torsion enumeration on tiny fields

## Stage F: Native order algorithms beyond exhaustive enumeration

### Why here

The current user-facing wrappers are good enough for the first milestone, but
they still rely on short-specific machinery for deeper algorithms in supported
characteristics.

### Deliverables

- separate which finite-field routes should remain shared/group-level from
  which ones should become model-owned
- decide whether quadratic-character counting belongs to the general model as a
  model-side cubic-in-`x` surface or remains a transported route
- decide separately for Schoof
- keep route-preserving reports explicit

### Dependencies

- Stage B required
- Stage E recommended before any native Schoof story

### Exit criteria

- the roadmap decision is encoded in code ownership, not only in prose
- any newly native route is tested against exhaustive counting on small fields

## Stage G: Graph layer and kernel witnesses owned by the general model

### Why here

Only after isogenies and torsion/kernel stories stabilize does it make sense to
let the graph layer stop treating the general model as a passenger.

### Deliverables

- define what a general-model kernel witness should store
- decide whether graph nodes store native general representatives, transported
  short representatives, or both
- add visualization/explanation support that says exactly which witness type is
  being shown

### Dependencies

- Stage D required
- Stage E strongly recommended

### Exit criteria

- at least one tiny graph workflow is owned by the general model in a way that
  is testable and inspectable

## Cross-cutting work

The following tasks should happen continuously, not only at the end:

### Documentation

- keep native-vs-transported boundaries explicit in rustdocs
- update module docs when a staged bridge becomes a native route
- keep this roadmap current when a stage is materially completed or reprioritized

### Visualization

- whenever a new public algebraic surface lands, add a matching educational
  formatter or explanation helper
- prefer side-by-side “general vs short companion” explanations when the bridge
  is still part of the story

### Testing

- preserve the dedicated `tests/compatibility.rs` contract
- add separate test files when one new stage grows its own vocabulary
- prefer exhaustive tiny-field checks before broad property tests when the
  represented regime allows it

### Examples

- keep one compact educational example per major public stage
- prefer examples that show what changes and what stays invariant under model
  transport

## Near-term recommendation

If work resumes immediately after this roadmap, the best next target is:

1. Stage A: native projective group law
2. Stage B: explicit admissible isomorphisms
3. Re-evaluate whether Stage C or Stage E should come next based on how much
   reuse the isomorphism layer actually unlocked

That sequence keeps the effort grounded in the two most reusable missing
foundations instead of jumping too early to Schoof or native isogenies.

## Non-goals for the next pass

The following should not be attempted all at once in one milestone:

- projective formulas
- function fields
- native Vélu
- division polynomials
- Schoof
- graph witnesses

That bundle is too large, too coupled, and too hard to validate in one pass.

