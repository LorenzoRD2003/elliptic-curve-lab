# AGENTS.md for `src/isogenies`

## Module mission

The `isogenies` module should introduce morphisms between elliptic curves
gradually, honestly, and in a way that helps the reader understand both the
mathematics and the implementation.

At this stage, the goal is not to ship a high-performance isogeny toolkit.
The goal is to make kernels, codomains, and point images easy to inspect and
easy to reason about in small finite examples.

## Educational posture

- Treat this module as educational infrastructure, not as production
  isogeny-based cryptography.
- Prefer explicit subgroup data, explicit curve data, and documented formulas
  over opaque black-box APIs.
- If an implementation only supports a classical setting, say so directly.
- If a formula is normalized, model-specific, or only valid in certain
  characteristics, document that fact where the API lives.
- If a piece of Vélu support is still partial, say so explicitly instead of
  pretending the full theory is already implemented.

## Design priorities

- Mathematical honesty before feature count.
- Small public APIs before broad generality.
- Explicit kernels before clever derived representations.
- Shared sources of truth before duplicated helper logic.
- Conservative, reviewable steps.

## Current posture

- Examples that exercise public isogeny workflows should require the
  `isogeny-lab` Cargo feature, and usually `visualization` as well when their
  main output is explanatory text. Keep this as a chapter-level feature rather
  than splitting immediately into `velu`, `graphs`, `frobenius`, or `dual`
  flags.
- The public `src/isogenies` root is now intentionally austere. Keep only
  genuinely reusable isogeny infrastructure there, and move
  short-Weierstrass-specific executable logic under
  `elliptic_curves::short_weierstrass::isogenies`.
- `IsogenyKernel` is a general finite-subgroup abstraction, not a Vélu-only
  helper.
- `KernelDescription` is the broader public kernel surface for general
  isogenies. Keep `IsogenyKernel` for the reduced explicit-subgroup case, and
  do not force inseparable kernels into a point-only representation.
- The current concrete construction path is Vélu for small finite subgroups on
  short-Weierstrass curves.
- Composition scaffolding is acceptable before full composed evaluation exists,
  as long as the docs say exactly what is implemented now and what still
  remains `todo!()`.
- `ComposedIsogeny` is now an actual small finite composition surface, not
  just a placeholder: strict composition, composition through an explicit
  bridge isomorphism, evaluation, degree, and exhaustive rational kernel
  materialization are all acceptable in the current implementation stage.
- Exhaustive equality helpers and exhaustive verification traits are now part
  of the intended confidence surface for this module. It is acceptable to
  compare maps by evaluating them on every rational point of a tiny curve, as
  long as the docs say so directly.
- Exhaustive dual search by enumerating tiny finite kernels and then checking
  the duality relations on rational points is acceptable for the current
  educational setting. Say directly that this is a small-field search
  routine, not a general dual-construction algorithm.
- The current short-Weierstrass dual search returns a concrete
  `DualVeluIsogeny<F>` represented as:
  - a Vélu isogeny on `E'`
  - followed by a base-field short-Weierstrass isomorphism back to `E`
- When those duality checks are really properties of one concrete
  `DualVeluIsogeny<F>` candidate relative to a chosen `φ`, prefer methods such
  as `dual.verify_left_dual_relation(phi)` and
  `dual.verify_right_dual_relation(phi)` over free helper functions.
- Generic Frobenius-vs-isogeny comparison reports such as
  `IsogenyFrobeniusRelation` and `IsogenyGraphFrobeniusReport` belong under
  `src/isogenies/frobenius_relation/`, not under
  `elliptic_curves::short_weierstrass::isogenies`, because they compare
  reusable isogeny/graph invariants rather than owning one concrete
  short-Weierstrass map construction.
- Within that `frobenius_relation/` subtree, prefer capability traits such as
  `FrobeniusComparableIsogeny` and `FrobeniusComparableIsogenyGraph` over
  free helper functions when the verification is naturally a property of the
  isogeny or graph object being checked.
- `ScalarMultiplicationIsogeny<C>` is an acceptable educational self-isogeny
  surface for small finite curves. Keep the docs explicit that
  `kernel_points()` means the rational kernel on `E(F_q)`, not the full
  geometric kernel over an algebraic closure.
- Scalar-multiplication isogenies, duality reports, Frobenius/Verschiebung
  degree summaries, and characteristic-factorization data should carry
  mathematical degrees/scalars as `BigUint`. Keep `usize` only for explicit
  enumerable kernel sizes or rational-point counts where the value is genuinely
  a collection length.
- Within `scalar_multiplication/`, prefer separating:
  - the isogeny type `[n]` itself
  - the characteristic-side factorization data `n = p^e m`
  - the kernel-structure analysis built from that factorization
  - the function-field / Verschiebung story
    rather than mixing those four narratives into one file.
- Frobenius-side certifications whose primary input is an explicit isogeny or
  an isogeny graph belong under `src/isogenies/frobenius_relation/`, even when the
  compared invariants are curve-order and trace data coming from
  `elliptic_curves::frobenius`.
- The current short-Weierstrass function-field pullback layer under
  `elliptic_curves::short_weierstrass::isogenies::function_field_maps` is an
  acceptable preparatory surface for later isogeny work. Represent a map
  `phi : E -> E'` by the pullback
  `phi^* : F(E') -> F(E)` through the images of `x'` and `y'`, and validate at
  least that those pullbacks live on `F(E)` and satisfy the codomain equation
  after substitution.
- The first separability-side surface should be a differential report rather
  than a bare boolean. It is acceptable for the current implementation to
  certify only the separable case from a non-zero differential multiplier and
  leave the zero-multiplier case as `Unknown` or another explicitly modest
  classification until the inseparable-factor machinery exists.
- That inseparable side now includes explicit absolute and relative Frobenius
  isogenies on short-Weierstrass curves:
  - absolute Frobenius lands on the `p`-power twist
  - relative Frobenius lands back on the same curve over `F_q`
  - both should expose their function-field pullbacks through the existing
    `ShortWeierstrassFunctionFieldMap` surface rather than through a parallel
    bespoke pullback type
  - both should report `separable_degree = 1` and a nontrivial
    `inseparable_degree`
  - their differential pullback reports should be reclassified as
    `PurelyInseparable`, not left in the generic zero-multiplier bucket
- A first `VerschiebungIsogeny` surface is acceptable even before full
  point-evaluation support, as long as the docs say clearly that it is
  presently a function-field-side witness certified by the composition
  relations with Frobenius and supplied `[p]` pullback maps, not yet an
  implementation of the full `Isogeny` trait.
- `ScalarMultiplicationIsogeny` on short-Weierstrass curves may now export its
  pullback map `[n]^*` directly by evaluating `[n]` on the generic point of
  `E(F(E))`.
- The certified Verschiebung route for `[p]^*` is still acceptable and useful
  as an independent characteristic-`p` certification path, but it is no longer
  the only available pullback construction for scalar multiplication.
- For curves over `F_{p^r}` that do not descend to `F_p`, reconstructing
  Verschiebung from `[p]^*` should be phrased as inverting the absolute
  Frobenius pullback `F(E^(p)) -> F(E)`, not as taking an ordinary `p`-th root
  inside one fixed function field.
- Keep the docs explicit that this pullback layer is presently weaker than a
  full certified isogeny constructor: it models the contravariant algebra map,
  but it does not yet prove that the data comes from a genuine finite isogeny
  or that the induced map is injective on function fields.
- For substitution into those pullback maps, prefer reusing the existing
  short-Weierstrass function-field arithmetic directly instead of introducing a
  second symbolic expression layer for `x` and `y`.
- Property-test fixtures for short-Weierstrass `function_field_maps` should
  generate genuinely valid pullback data, not arbitrary pairs of functions.
  The current acceptable families are:
  - self-maps such as identity or `y -> -y` on one curve
  - constant maps to rational finite codomain points
  - composable chains built from those same valid map families
- Exhaustive map comparison helpers such as pointwise equality on `E(F_q)` are
  acceptable under `src/isogenies/` when they improve confidence in small
  examples, dual checks, or composition checks.
- Kernel validation is intentionally exhaustive for small groups:
  identity, on-curve membership, closure under negation, and closure under
  addition.
- The current `VeluIsogeny<C>` design uses the same Rust model type for domain
  and codomain as an implementation simplification. This means “same curve model
  family”, not “same curve value”.
- The `graphs/` subtree is now the educational `ℓ`-isogeny graph surface for
  small `Fp<P>` curves. Keep it explicit, representative-based, and gradual
  instead of letting it sprawl into a large generic graph framework.
- Graph edges may store explicit target-representative transport, but prefer a
  model-level associated isomorphism witness over field-specific plumbing so
  the representation can later grow beyond one concrete curve family.
- Within `graphs/`, prefer keeping the public barrel type-first and namespace
  austere: reexport the main graph/report/value types from `graphs/mod.rs`,
  and feel free to lower helper submodules such as `edge`, `node`,
  `verification`, or `builder` to `pub(crate)` once their public types are
  already reexported at the top level.
- While the graph container stays vector-backed, keep the node-id contract
  simple and explicit: dense ids may follow vector order, and graph-side lookup
  helpers should stay honest about that assumption.
- As the graph subtree grows, it is acceptable to split narrowly focused
  helper logic into files such as `torsion.rs` instead of overloading
  `builder.rs` with unrelated responsibilities.
- However, keep the ownership boundary explicit:
  - generic exact-order point logic belongs under `elliptic_curves`
  - `graphs/torsion.rs` should now be a thin wrapper specialized to graph-side
    cyclic-kernel construction, ideally exposed through one narrow graph-side
    capability trait on suitable curve models rather than through a public
    free function
  - division-polynomial torsion recovery may feed graph/kernel
    workflows, but the recovery logic itself still belongs under
    `elliptic_curves::short_weierstrass::division_polynomials`
- Now that the graph subtree has both container logic and construction logic,
  keep `graphs/builder/` as a small module tree, for example separating
  generic graph storage, short-Weierstrass-specific construction, and
  `builder/tests.rs`.
- For graph-side heuristic reports such as `VolcanoLikeLayering`, prefer small
  semantic helpers like `role_of`, `nodes_at_level`, or `count_role` over
  making every consumer re-scan raw `(node, role)` slices by hand.
- Likewise, prefer keeping `IsogenyGraphNode` and `IsogenyGraphEdge` public as
  lightweight structural carriers (`id`, endpoints, degree, counts), while
  lowering storage-heavy details such as representative curves or explicit
  kernel subgroups to `pub(crate)` once visualization and internal algorithms
  can consume them inside the crate.
- After that split, keep the generic graph container and traversal helpers
  under `isogenies::graphs`, but keep short-Weierstrass graph-construction
  implementations attached to the graph/builder owner that actually uses
  them, instead of creating a dedicated `graphs::short_weierstrass` namespace
  unless that subtree grows a second genuinely independent story.
- For graph structural analysis, helpers such as weakly connected
  components belong under `src/isogenies/graphs/` rather than under
  `visualization/`. Presentation layers should reuse those helpers instead of
  recomputing graph structure independently.
- The same rule applies to directed-cycle helpers: keep cycle discovery under
  `src/isogenies/graphs/` as structural analysis, and let visualization layers
  consume those helpers instead of embedding graph traversal logic of their own.
- For the first graph-side endomorphism report surface, prefer one aggregate
  public entry point such as `IsogenyGraph::endomorphism_report_at(&ell)`.
  Keep local edge classifiers, candidate-set comparison helpers, and
  graph-vs-volcano comparison reports internal until an example or
  visualization needs them as part of a stable educational story.
- In graph-side aggregate reports, it is acceptable for per-edge report values
  to repeat structural identifiers such as `source` and `target` even when
  those identifiers also live on the underlying graph edge. That duplication
  keeps the report independently readable and avoids forcing educational
  consumers to bounce back to the graph just to understand one edge annotation.
- Generic orchestration belongs in `velu/core.rs`; model-specific formulas
  belong in the specialized short-Weierstrass subtree under
  `elliptic_curves::short_weierstrass::isogenies::velu`.
- Internal precomputations such as `VeluKernelData` should remain the shared
  source of truth for codomain and evaluation logic when they arise from the
  same kernel formulas.
- The current short-Weierstrass Vélu layer now also exports the rational map on
  function fields. Keep those pullbacks synchronized with the same
  translation-sum normalization already used for point evaluation:
  `x_phi(P) = x(P) + sum_Q (x(P + Q) - x(Q))` and
  `y_phi(P) = y(P) + sum_Q (y(P + Q) - y(Q))`.

## Kernel rules

- A kernel is a finite subgroup of the domain curve, not merely an arbitrary
  list of points.
- Public constructors should keep that meaning explicit:
  - `new(...)` validates an explicit finite subgroup
  - `cyclic(...)` enumerates the subgroup generated by one point in a small
    finite setting
- Prefer `HashSet`-based input when the mathematical object is a set rather
  than a sequence.
- Preserve and document invariants that downstream code relies on. The current
  example is that the identity is stored first in `IsogenyKernel::points()`.
- Do not move kernel logic into `velu/` unless the change reflects a real
  mathematical restriction rather than just local convenience.

## Vélu guidance

- Keep the conceptual separation clear between:
  - the abstract idea of an isogeny with kernel `G`
  - the specific Vélu construction from a finite subgroup `G`
- Prefer constructors that start from the mathematical input a learner
  naturally has, such as:
  - a generator of a cyclic subgroup
  - an explicit finite subgroup
- Generic constructors like `from_points(...)`, `from_generator(...)`, and
  `from_kernel(...)` should stay in the generic core when that remains honest.
- Model-specific codomain or evaluation formulas should not leak back into the
  generic core.
- Model-specific dual search for Vélu currently belongs under
  `elliptic_curves::short_weierstrass::isogenies::velu::dual` rather than the
  generic core, since it depends on short-Weierstrass isomorphism witnesses and
  small-field exhaustive search.
- If codomain formulas and evaluation formulas come from the same
  normalization, keep them coupled through shared internal data instead of
  duplicating derivations in two places.
- If a helper computes affine coordinates and returns `None`, document whether
  that means “no affine image because the image is the identity” rather than
  “the computation failed”.

## Scope guidance

- Small finite examples are the intended proving ground for this module.
- Exhaustive or direct algorithms are acceptable when the docs say so clearly.
- It is fine to support short-Weierstrass curves first.
- For graph work, keep the initial scope explicitly limited to:
  - short-Weierstrass curves over small `Fp<P>`
  - rational cyclic kernels generated by points of order `ℓ` in `E(F_p)`
  - codomain deduplication up to base-field isomorphism, not algebraic-closure
    isomorphism
  - volcanoes as educational heuristics, not a complete theory of quadratic
    orders
- Even while the first implementation is short-Weierstrass-first, it is
  acceptable to shape the graph data around generic curve-model capabilities
  such as a `j`-invariant trait plus a same-family isomorphism witness type.
- For node deduplication in graph, it is acceptable to use a
  short-Weierstrass-specific helper that filters first by equal `j` and then
  confirms base-field isomorphism by exhaustive witness search on the small
  field.
- For exact-order point checks in graph, it is acceptable and often
  clearer to use the explicit “`[n]P = O` plus no `[(n/p)]P = O` for prime
  divisors `p | n`” criterion directly.
- Now that that exact-order logic exists generically under `elliptic_curves`,
  graph code should reuse it instead of reimplementing it locally.
- For kernel and torsion enumeration in the graph layer, deduplicate cyclic kernels
  by subgroup equality, not by chosen generator. In particular, `P` and `-P`
  should collapse to the same cyclic kernel whenever they generate the same
  subgroup.
- For graph outgoing-edge construction, prefer the pipeline
  “exact-order rational points -> distinct cyclic kernels -> Vélu codomains ->
  node deduplication by base-field isomorphism” instead of deriving edges from
  raw torsion generators directly.
- It is acceptable for division-polynomial tests to verify that
  division-polynomial-derived exact torsion generators feed this same cyclic
  kernel and Vélu pipeline without duplicating kernel logic locally.
- When that outgoing-edge pipeline is implemented, it is reasonable to factor
  the “exact target / isomorphic target / insert fresh target” decision into a
  small private resolver helper rather than repeating it inline inside the
  edge-building loop.
- When a newly produced Vélu codomain already matches an existing stored
  representative exactly, it is acceptable and clearer to record
  `EdgeTargetWitness::Identity`; reserve explicit witnesses for the genuine
  “same isomorphism class, different stored representative” case.
- For whole-graph construction in this layer, breadth-first expansion by a
  bounded edge depth is a good educational default. Treat nodes without
  rational cyclic kernels of the requested prime degree as ordinary leaves, not
  as a global build failure.
- For local verification in this layer, it is reasonable to reconstruct each
  stored edge as “Vélu from the stored kernel, followed by the stored target
  transport witness” and then reuse the dual-isogeny exhaustive checks on rational
  points. Prefer checking reverse-direction edges already present in the graph
  before considering any fresh global dual search.
- Current technical debt: the first educational graph summaries treat the graph
  as uniform-degree and read one `ℓ` value from the stored edge set. If a later
  phase allows mixed-degree edges in one graph container, the summary surface
  must be widened explicitly instead of silently reusing one scalar `degree`
  field.
- For volcano helpers in the graph layer, keep the first implementation explicitly
  graph-theoretic and root-dependent: weak BFS layering plus local weak degree
  is acceptable as an educational heuristic, but the docs should say directly
  that this does not compute endomorphism rings, certify ordinary components,
  or prove a true Kohel/Sutherland volcano structure.
- If a later bridge report compares that heuristic with arithmetic data from
  `elliptic_curves::endomorphisms`, keep the ownership split explicit:
  arithmetic candidate orders and local levels still belong to
  `elliptic_curves`, while the comparison report lives under
  `isogenies::graphs` as a consumer of both layers.
- Such a bridge report should stay honest and modest: comparing candidate
  local levels with heuristic BFS layers is acceptable, but do not present the
  result as a certification of ascending, descending, or horizontal edges, nor
  as a proved volcanic level of the curve.
- The same caution applies to tentative edge-level endomorphism reports:
  names like `PossiblyHorizontal`, `PossiblyAscending`, and
  `PossiblyDescending` are acceptable when they arise only from comparing
  candidate local levels, but keep `Ambiguous` and `Unsupported` available and
  say directly that the result is not a definitive edge classification.
- If graph nodes already store concrete finite-field curve representatives,
  it is acceptable to derive each node's `EndomorphismRingCandidateSet`
  automatically from the stored curve via the existing Frobenius pipeline,
  rather than asking callers to provide those candidate sets manually.
- If the graph layer grows several endomorphism-side bridge files, prefer
  grouping them under a dedicated internal submodule such as
  `isogenies::graphs::endomorphisms` rather than leaving them flat beside core
  graph storage files like `node.rs`, `edge.rs`, and `volcano.rs`.
- If a later aggregate graph-side endomorphism report is added, prefer
  building it from:
  `node -> candidate set -> local levels`
  together with the already implemented tentative edge relations, rather than
  introducing a second independent source of truth for node arithmetic.
- Prefer names that keep the heuristic status visible in the API surface. A
  name such as `VolcanoLikeLayering` is better than something that sounds like
  a canonical arithmetic decomposition.
- Do not let the volcano-like heuristic grow without a clear boundary between
  visual intuition and mathematical claim. If a future extension starts to
  sound like certification rather than pedagogy, either strengthen the
  mathematics explicitly or rename/re-scope the surface so the educational
  status remains obvious.
- A runnable graph example may present the whole graph workflow in one
  short script: build the graph, print a structural summary, run local
  verification, and then explain one chosen volcano-like layering. Prefer one
  fixed small-field curve over an internal search so the example stays short
  and reproducible.
- Test graph behavior at the graph level, not only through summaries: for
  example, it is worth checking explicitly that distinct rational kernels may
  survive as parallel edges to one deduplicated codomain node, that reverse
  edges appear in known small examples, and that every stored edge can be
  reconstructed back into a valid local Vélu-based map.
- Do not let graph tests rely only on `ℓ = 2`; keep at least one small graph-side
  example with `ℓ = 3` or `ℓ = 5` so the graph layer exercises non-binary
  behavior too.
- When graph-side torsion wrappers coexist with division-polynomial
  tooling, include at least one agreement test showing that graph cyclic
  kernels of order `ℓ` match the exact-order torsion recovered from division
  polynomials on a small sample curve.
- If the heuristic needs an “I cannot classify this from the chosen root”
  escape hatch, prefer an explicit `Unknown` role over pretending every stored
  node belongs to one inferred volcano level.
- Use the term `reverse edge` for graph adjacency facts and reserve
  `dual` for the stronger verified property coming from the two duality
  relations.
- When the rational `ℓ`-torsion splits as `E[ℓ](F_q) ≅ Z/ℓZ × Z/ℓZ`, the
  educational graph layer should be prepared to observe `ℓ + 1` distinct
  cyclic kernels of order `ℓ`.
- It is also fine for the public `Isogeny` trait to stay minimal while the
  internal Vélu support grows underneath it.
- For composed isogenies, prefer field names like `first` and `second` over
  `left` and `right` when that reduces right-to-left composition confusion for
  readers.
- Do not jump early to large-degree optimizations, kernel polynomials,
  square-root Vélu, modular-polynomial navigation, or cryptographic hardening.
- Do not introduce new abstraction layers unless they remove real duplication
  or mark a real mathematical boundary.

## Error conventions

- Keep recoverable isogeny-domain failures in `IsogenyError`.
- Prefer keeping `IsogenyError` as a small top-level wrapper around narrower
  sub-enums such as kernel, verification, map/function-field, duality, and
  Verschiebung errors, instead of growing one flat catch-all enum again.
- Reuse `CurveError` through `.into()` `IsogenyError::Curve(...)` when the failure is
  really a curve-domain issue.
- Prefer specific subgroup-validation failures over generic “invalid kernel”
  errors when the distinction teaches something real.
- If a characteristic restriction matters for the chosen formulas, expose it as
  a typed error rather than a string.
- The graph subtree may introduce its own typed `IsogenyGraphError` surface
  when graph-building failures are distinct from lower-level curve or isogeny
  failures. Prefer wrapping `CurveError` and `IsogenyError` explicitly instead
  of collapsing them into strings.

## Testing expectations

- Test both valid and invalid kernel construction.
- Test cyclic kernels on at least one exact small example.
- Test the identity case for evaluation explicitly.
- Test that kernel points map to the codomain identity.
- Test at least one non-kernel point image on a small exact example.
- When a codomain formula is implemented, test that the resulting curve is the
  expected one on a concrete small example.
- For composition scaffolds, test both the accepted exact-middle-curve case and
  the rejected mismatched-middle-curve case before adding richer behavior.
- When evaluation and codomain construction are both implemented, add tests for
  their coherence rather than testing them in isolation only.
- As the module matures, add structural isogeny tests such as:
  - constancy on kernel cosets
  - homomorphism behavior
  - cardinality relations in small finite settings
- For scalar-multiplication isogenies, test the current `n != 0` policy
  explicitly. If `[1]` is used as an identity-like map in composition tests,
  make that role explicit in the test name and docs.
- For dual search, prefer at least one degree-2 and one degree-3 example over
  small prime fields, and test both the existence of the dual and the expected
  rational-kernel size of the returned candidate.
- For graph scaffolding, prefer at least one medium small-field
  test such as the existing `F41` Vélu example instead of relying only on
  very tiny fields like `F7`.
- Test typed error variants directly.

## Documentation expectations

- Public items should explain the mathematical object they represent.
- If a construction assumes characteristic different from `2` and `3`, say so
  directly.
- If a helper is only honest for small finite groups, say so directly.
- If the codomain is computed by a specific normalized version of Vélu's
  formulas, document the exact formulas being used.
- If a composed or dual isogeny depends on a bridge isomorphism, say directly
  whether the bridge is strict identity, an explicit cached witness, or an
  exhaustively chosen witness among several valid base-field isomorphisms.
- If a dual isogeny is represented as “Vélu part plus isomorphism back to the
  original curve”, say so directly rather than describing it as an opaque
  black-box dual constructor.
- If the current implementation makes a capability-level simplification, such as
  keeping domain and codomain in the same Rust model family, document that as a
  simplification rather than as a universal theorem.
- Use concrete subgroup notation such as `G = <P>` and concrete examples when
  they make the implementation easier to follow.

## Review heuristics

A good change under `src/isogenies` should improve at least one of:

- kernel invariants
- formula coherence
- readability
- mathematical honesty
- test coverage

If a change makes it harder to explain what the kernel is, how the codomain is
obtained, or why an image point should land on that codomain, it is probably
moving too fast for the current phase of the project.
