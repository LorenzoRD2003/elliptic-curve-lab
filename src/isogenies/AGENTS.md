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

- `IsogenyKernel` is a general finite-subgroup abstraction, not a Vélu-only
  helper.
- The current concrete construction path is Vélu for small finite subgroups on
  short-Weierstrass curves.
- Composition scaffolding is acceptable before full composed evaluation exists,
  as long as the docs say exactly what is implemented now and what still
  remains `todo!()`.
- `ComposedIsogeny` is now an actual small finite composition surface, not
  just a placeholder: strict composition, composition through an explicit
  bridge isomorphism, evaluation, degree, and exhaustive rational kernel
  materialization are all acceptable in the current milestone.
- Exhaustive equality helpers and exhaustive verification traits are now part
  of the intended confidence surface for this module. It is acceptable to
  compare maps by evaluating them on every rational point of a tiny curve, as
  long as the docs say so directly.
- Exhaustive dual search by enumerating tiny finite kernels and then checking
  the duality relations on rational points is acceptable for the current
  educational milestone. Say directly that this is a small-field search
  routine, not a general dual-construction algorithm.
- The current short-Weierstrass dual search returns a concrete
  `DualVeluIsogeny<F>` represented as:
  - a Vélu isogeny on `E'`
  - followed by a base-field short-Weierstrass isomorphism back to `E`
- Public helpers such as `verify_left_dual_relation(...)` and
  `verify_right_dual_relation(...)` are acceptable when they make the dual
  identities explicit and reusable in tests, examples, and visualization.
- `ScalarMultiplicationIsogeny<C>` is an acceptable educational self-isogeny
  surface for small finite curves. Keep the docs explicit that
  `kernel_points()` means the rational kernel on `E(F_q)`, not the full
  geometric kernel over an algebraic closure.
- Exhaustive map comparison helpers such as pointwise equality on `E(F_q)` are
  acceptable under `src/isogenies/` when they improve confidence in small
  examples, dual checks, or composition checks.
- Kernel validation is intentionally exhaustive for small groups:
  identity, on-curve membership, closure under negation, and closure under
  addition.
- The current `VeluIsogeny<C>` design uses the same Rust model type for domain
  and codomain as a milestone simplification. This means “same curve model
  family”, not “same curve value”.
- The upcoming `graphs/` subtree is intended for milestone-6 educational
  `ℓ`-isogeny graph exploration over small `Fp<P>` curves. Start with a small
  scaffold and fill it gradually instead of jumping directly to a large graph
  framework.
- Graph edges may store explicit target-representative transport, but prefer a
  model-level associated isomorphism witness over field-specific plumbing so
  the representation can later grow beyond one concrete curve family.
- While the graph container stays vector-backed, keep the node-id contract
  simple and explicit: dense ids may follow vector order, and graph-side lookup
  helpers should stay honest about that assumption.
- As the graph subtree grows, it is acceptable to split narrowly focused
  helper logic into files such as `torsion.rs` instead of overloading
  `builder.rs` with unrelated responsibilities.
- Once the graph subtree has both container logic and construction logic, it
  is acceptable and preferable to split `graphs/builder.rs` into a small
  `graphs/builder/` module tree, for example separating generic graph storage,
  short-Weierstrass-specific construction, and `builder/tests.rs`.
- For milestone-6 structural graph analysis, helpers such as weakly connected
  components belong under `src/isogenies/graphs/` rather than under
  `visualization/`. Presentation layers should reuse those helpers instead of
  recomputing graph structure independently.
- The same rule applies to directed-cycle helpers: keep cycle discovery under
  `src/isogenies/graphs/` as structural analysis, and let visualization layers
  consume those helpers instead of embedding graph traversal logic of their own.
- Generic orchestration belongs in `velu/core.rs`; model-specific formulas
  belong behind the private support trait and the specialized subtree under
  `velu/short_weierstrass/`.
- Internal precomputations such as `VeluKernelData` should remain the shared
  source of truth for codomain and evaluation logic when they arise from the
  same kernel formulas.

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
  `velu/short_weierstrass/` rather than the generic core, since it depends on
  short-Weierstrass isomorphism witnesses and small-field exhaustive search.
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
- For milestone 6 graph work, keep the initial scope explicitly limited to:
  - short-Weierstrass curves over small `Fp<P>`
  - rational cyclic kernels generated by points of order `ℓ` in `E(F_p)`
  - codomain deduplication up to base-field isomorphism, not algebraic-closure
    isomorphism
  - volcanoes as educational heuristics, not a complete theory of quadratic
    orders
- Even while the first implementation is short-Weierstrass-first, it is
  acceptable to shape the graph data around generic curve-model capabilities
  such as a `j`-invariant trait plus a same-family isomorphism witness type.
- For node deduplication in milestone 6, it is acceptable to use a
  short-Weierstrass-specific helper that filters first by equal `j` and then
  confirms base-field isomorphism by exhaustive witness search on the small
  field.
- For exact-order point checks in milestone 6, it is acceptable and often
  clearer to use the explicit “`[n]P = O` plus no `[(n/p)]P = O` for prime
  divisors `p | n`” criterion directly.
- For milestone-6 kernel and torsion enumeration, deduplicate cyclic kernels
  by subgroup equality, not by chosen generator. In particular, `P` and `-P`
  should collapse to the same cyclic kernel whenever they generate the same
  subgroup.
- For milestone-6 outgoing-edge construction, prefer the pipeline
  “exact-order rational points -> distinct cyclic kernels -> Vélu codomains ->
  node deduplication by base-field isomorphism” instead of deriving edges from
  raw torsion generators directly.
- When that outgoing-edge pipeline is implemented, it is reasonable to factor
  the “exact target / isomorphic target / insert fresh target” decision into a
  small private resolver helper rather than repeating it inline inside the
  edge-building loop.
- When a newly produced Vélu codomain already matches an existing stored
  representative exactly, it is acceptable and clearer to record
  `EdgeTargetWitness::Identity`; reserve explicit witnesses for the genuine
  “same isomorphism class, different stored representative” case.
- For milestone-6 whole-graph construction, breadth-first expansion by a
  bounded edge depth is a good educational default. Treat nodes without
  rational cyclic kernels of the requested prime degree as ordinary leaves, not
  as a global build failure.
- For milestone-6 local verification, it is reasonable to reconstruct each
  stored edge as “Vélu from the stored kernel, followed by the stored target
  transport witness” and then reuse the M5 exhaustive checks on rational
  points. Prefer checking reverse-direction edges already present in the graph
  before considering any fresh global dual search.
- Current technical debt: the first educational graph summaries treat the graph
  as uniform-degree and read one `ℓ` value from the stored edge set. If a later
  phase allows mixed-degree edges in one graph container, the summary surface
  must be widened explicitly instead of silently reusing one scalar `degree`
  field.
- For milestone-6 volcano helpers, keep the first implementation explicitly
  graph-theoretic and root-dependent: weak BFS layering plus local weak degree
  is acceptable as an educational heuristic, but the docs should say directly
  that this does not compute endomorphism rings, certify ordinary components,
  or prove a true Kohel/Sutherland volcano structure.
- Do not let the volcano-like heuristic grow without a clear boundary between
  visual intuition and mathematical claim. If a future extension starts to
  sound like certification rather than pedagogy, either strengthen the
  mathematics explicitly or rename/re-scope the surface so the educational
  status remains obvious.
- A milestone-6 runnable example may present the whole graph workflow in one
  short script: build the graph, print a structural summary, run local
  verification, and then explain one chosen volcano-like layering. Prefer one
  fixed small-field curve over an internal search so the example stays short
  and reproducible.
- Test graph behavior at the graph level, not only through summaries: for
  example, it is worth checking explicitly that distinct rational kernels may
  survive as parallel edges to one deduplicated codomain node, that reverse
  edges appear in known small examples, and that every stored edge can be
  reconstructed back into a valid local Vélu-based map.
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
- Reuse `CurveError` through `IsogenyError::Curve(...)` when the failure is
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
- For milestone-6 graph scaffolding, prefer at least one medium small-field
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
- If the current implementation makes a milestone-level simplification, such as
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
