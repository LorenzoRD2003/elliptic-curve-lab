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
- The public `crate::visualization` module is compiled behind the
  `visualization` Cargo feature. Examples that depend on this subtree should
  declare `required-features = ["visualization"]`, plus any chapter feature
  such as `analytic`, `isogeny-lab`, or `advanced-point-counting` that
  describes the mathematical story being shown.
- Domain-layer accessors that exist only so visualization can explain an
  internal report should be gated with `#[cfg(feature = "visualization")]`
  rather than compiled unconditionally and left as dead code in non-visual
  builds.

## Structure rules

- `visualization/fields/` is for field-domain values and field-domain
  explanations.
  That includes quotient-style field values such as polynomial quotients and
  rational functions, plus short summaries of their ambient field families.
- `visualization/elliptic_curves/` is for curve equations, points, group-law
  explanations, small finite curve-group reports, and division-polynomial /
  rational-torsion summaries.
  It also hosts analytic summaries for lattices, Eisenstein
  truncations, analytic invariants, torus-to-curve maps, and differential
  equation checks.
- `visualization/isogenies/` is for educational summaries of kernels,
  codomains, point-evaluation formulas for explicit isogeny constructions, and
  summaries of composition, scalar multiplication, dual
  isogenies, exhaustive dual verification on tiny curves, and
  graph summaries and adjacency explanations.
  Prefer file names that say which family of summaries they hold, such as
  `velu`, `graph`, or `derived_maps`, over vague catch-all names.
- `visualization/polynomials/` is for polynomial-domain values and
  polynomial-domain explanations.
- If a helper explains a capability trait such as `SqrtField`, it belongs in
  the matching mathematical subtree rather than in a generic catch-all file.
- If a new algorithmic route introduces its own report type in the domain
  layer, prefer adding a dedicated visualization helper for that route-level
  report instead of explaining it only indirectly through one larger unified
  enum wrapper.
- When a runnable example prints one report or value that implements
  `Visualizable`, prefer calling `.describe()` at the example call site.
  Reserve route-specific `explain_*` helpers for lists, context-dependent
  explanations, or cases where no single stored value owns the whole story.
- Keep the public root API of `crate::visualization` narrow. Prefer exporting
  visualization traits such as `Visualizable` over reexporting route-specific
  `describe_*`, `explain_*`, or `format_*` helpers. Helpers may remain inside
  the private visualization implementation tree while examples and callers use
  `.describe()` / `.format_compact()` on value objects.
- Helpers that back `Visualizable` implementations should not be public API.
  Keep them private when used in one file, or use `pub(crate)` when sibling
  visualization modules need to share formatting glue.
- If the same small formatting helper appears in several visualization files,
  move it to `visualization::shared` instead of cloning local definitions.
- For plain Boolean fields rendered as English yes/no text, use
  `visualization::shared::yes_no` instead of repeating local `if ... { "yes" }
  else { "no" }` expressions.
- For comma-separated lists of already formatted values, prefer
  `visualization::shared::comma_list`; for comma-separated compact
  `Visualizable` values, prefer `compact_visualizable_list`.
- When one visualization file grows several independent report families, split
  it into private sibling submodules by responsibility. Keep only the narrow
  `pub(crate)` reexports needed by other visualization modules, and leave
  caller-facing examples on the `.describe()` surface.

## Honesty rules

- Say when arithmetic is exact.
- Say when arithmetic is approximate.
- Say when a report depends on one or more truncation radii and which
  quantities were approximated through those truncations.
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
- If the domain API already separates disjoint cases through a typed enum,
  prefer mirroring that split in visualization helpers instead of flattening
  everything back into one generic textual shape.

## Formatting guidance

- Prefer deterministic plain text.
- Prefer English wording in user-facing output unless a task explicitly asks
  for another language.
- For compact numeric formatting, prefer dropping pointless trailing zeros and
  simplifying special coefficients such as `1 + 1i` into `1 + i` when that
  improves readability without hiding meaningful scale.
- For rational-function helpers, prefer reusing the existing dense-polynomial
  formatters for numerators and denominators so quotient notation stays
  visually consistent with the polynomial layer.
- Tables are welcome when they genuinely help, especially for small finite
  fields.
- Explanations should highlight the important intermediate quantities, not
  every possible low-level detail.
- For isogenies, prefer showing the kernel points, the codomain formulas, and
  a few key translation-sum terms over dumping large algebraic expressions
  without guidance.
- For labeled crater-walk summaries, keep graph diagnostics and arithmetic
  labels visibly separate: show the local ideal/form labels and the
  graph-deterministic walk, but explicitly say that arithmetic orientation and
  class-group action are not certified yet.
- For crater reports, keep the basic certified-evidence summary on
  `CraterReport::describe()` instead of hand-formatting crater prime, node
  count, shape, and horizontal-cycle count in examples.
- For short-Weierstrass function-field pullbacks `phi^* : F(E') -> F(E)`,
  prefer showing:
  - the direction `F(E') -> F(E)` explicitly
  - the stored generator images `phi^*(x')` and `phi^*(y')`
  - the substitution rule for `A(x')`, `r(x')`, or `A(x') + y'B(x')`
  - the contravariant composition rule `(psi o phi)^* = phi^* o psi^*`
- For differential pullback reports attached to those maps, prefer showing:
  - `ω_E = dx/(2y)` and `ω_E' = dx'/(2y')`
  - the intermediate derivative `dX_φ/dx`
  - the factor of `dx` in `φ^*(ω_E')`
  - the multiplier `c_φ = y*(dX_φ/dx)/Y_φ`
  - the current separability classification in plain language
- For graph visualizations, prefer one compact structural summary plus explicit
  node/edge listings and adjacency lists. Say directly that nodes store
  representatives and edges may carry transport witnesses onto those stored
  representatives.
- For graph-side endomorphism candidate refinement visualizations, say that
  surviving orders are compatible with observed evidence at `ℓ`; do not present
  a unique survivor as a certificate of the exact endomorphism ring `End(E)`.
- For graph-side endomorphism-ring level recovery visualizations, keep local
  and global stories separate: local reports explain `e`, `δ`, and
  `d = e - δ` at one prime `ℓ`, while global reports only assemble already
  supplied local reports into `u = ∏ℓ^{d_ℓ}` when every prime divisor of `v`
  is covered. If evidence is partial, show the missing primes instead of
  implying that `End(E)` has been fully recovered. When a runnable example
  demonstrates global assembly and can do so cheaply, prefer showing both a
  partial assembly and the completed multi-prime assembly. Keep visualization
  helpers presentation-only over existing reports; graph construction and
  multi-prime recovery orchestration should go through the graph-side API such
  as `recover_endomorphism_ring_at(...)`.
- For class-group-action intro visualizations, keep horizontal ideal reports
  framed as compatibility between certified crater evidence and a supplied
  prime-norm ideal. Crater walk reports should be framed as deterministic
  graph walks through certified crater evidence, not as oriented class-group
  actions. Do not describe either surface as a computation of `[𝔞] * E`.
- For analytic output, prefer showing:
  - the chosen `τ` or lattice basis
  - the derived modular parameter `q = e^{2π i τ}` when a routine is expressed
    through `q`-expansions
  - the modular matrix or accumulated modular matrix when modular actions or
    fundamental-domain reductions are being explained
  - the truncation radii
  - the approximate complex values actually computed
  - whether a comparison held approximately, failed, or hit a pole
  - when two analytic routes are compared, both approximations and the
    residual `difference` / `|difference|`
- For current analytic period-recovery output, prefer showing:
  - the chosen or recovered period basis `ω₁, ω₂` and the implied modulus `τ`
  - when a canonical modular representative is available, both the natural
    recovered `τ` and the canonically reduced one
  - the recovered cubic roots in their stored order, while saying that the
    order is not canonical
  - the reconstructed invariants or `j`-comparison residuals
  - for inverse-uniformization checks, the validation `τ`, the recovered
    invariant side, the curve-side `j`, and the truncation radius used to
    recompute the analytic lattice invariants
  - when scale-sensitive invariant checks are present, whether the outcome was
    direct agreement, only same modular class via `j`, or a genuinely
    inconsistent recovery
  - the numerical method/status and the phase counters that explain where the
    work went
  - any Cardano-branch diagnostics that explain why one branch pair was chosen
  - any geometric root classification separately from any “nearly repeated”
    warning
  - for Legendre reduction, the chosen `λ`, the six-element orbit when useful,
    the selected permutation, and the distinction between the branch-independent
    right-hand-side scale factor and the principal-branch `y` / differential scales
  - for canonicalized `τ`, the accumulated modular matrix and whether the
    reduction was already reduced or required actual modular steps
  - when `λ` is near `0`, `1`, or `∞`, or when distinct roots would collapse
    to the same compact string, prefer a higher-precision diagnostic formatter
    over the compact pretty-printer
- For division-polynomial explanations, prefer showing:
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
- For dual-isogeny summaries, it is acceptable to use compact symbolic
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
- For division-polynomial helpers, test both the compact summary surface and the
  richer torsion explanation surface for the phrases the runnable example
  relies on.
- For analytic helpers, test both the finite and pole-facing branches when
  applicable, and verify the phrases the runnable analytic examples rely on.
- When a helper is intended to support a runnable example, test the important
  phrases that the example relies on rather than snapshotting the entire final
  console output.
- For staged model-introduction milestones, pair each new visualization helper
  with a runnable example that shows the model itself, any explicit reduction
  witness available at this stage, and one side-by-side computation that
  confirms the educational story.

## Review heuristics

A good change under `src/visualization` should improve at least one of:

- readability
- mathematical honesty
- pedagogical usefulness
- consistency with the underlying domain API

If a visualization helper hides an important caveat from the reader, it is not
done yet.
