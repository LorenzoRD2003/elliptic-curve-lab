# AGENTS.md for `src/numerics`

## Module mission

The `numerics` module hosts small shared numerical infrastructure that is more
specialized than `utils` but broader than any single mathematical domain.

Its job is to keep approximate policies explicit and reusable across modules
such as `fields` and `elliptic_curves`.

## Design priorities

- Prefer tiny, well-documented value objects.
- Make numerical policy explicit rather than hidden in global constants.
- Keep APIs educational and honest about approximation.
- Avoid turning this module into a grab bag of unrelated helpers.

## Scope guidance

- Good fits include tolerances, truncation settings, normalization choices,
  and similarly small numerical policy types.
- Reusable geometry for numerical work in the complex plane, such as line
  segments, rays, or compactifying parameterizations used by several analytic
  routines, may also live here when they do not depend on one specific
  integrand or elliptic-curve convention.
- Generic quadrature helpers such as composite Simpson integration also
  belong here when they depend only on interval data and sampled values, and
  do not encode one domain-specific branch-choice or contour convention.
- If the even-budget normalization is only an implementation detail of the
  quadrature domain object, prefer keeping it internal to the module instead
  of exposing a second standalone public helper for the same policy.
- When a public numerical error enum is reached only through already-validated
  value objects, avoid keeping variants for impossible pre-validation states.
  Public error surfaces should reflect reachable failures from the exposed API.
- When a shared quadrature helper repeatedly takes the same interval endpoints
  plus a requested budget, prefer a small validated value object for that
  domain instead of passing three loose scalars through every call surface.
- If both a domain-based API and scalar-based convenience wrappers exist
  temporarily, prefer converging back to the domain-based surface once the
  value object proves useful, so the long-term public API stays smaller.
- When a shared quadrature rule has a classical weight pattern, prefer
  factoring that pattern into a named helper instead of leaving it buried in
  one loop-local conditional. That makes the formula easier to audit against
  the mathematics.
- When such a helper only computes one node's weighted contribution, prefer a
  pure return-value helper over a side-effecting “mutate the running sum”
  helper unless mutation is genuinely the clearer interface.
- When such complex-path geometry is exposed here, keep it integrand-agnostic
  and deterministic: endpoints, directions, affine interpolation, and
  compactifying parameterizations are good fits; branch-choice policy,
  singular-locus avoidance, or Abel-Jacobi-specific heuristics are not.
- If a compactified path parameterization is exposed here, also prefer
  exposing its derivative when downstream numerical integration would
  otherwise have to duplicate that calculus formula in domain-specific code.
- Likewise, if downstream code needs the inverse of that compactification,
  prefer exposing the inverse map here as part of the same geometric surface
  instead of recomputing the rational formula ad hoc elsewhere.
- Small exact numerical sequences or coefficient helpers that are shared by
  several mathematical domains may also live here, as long as they are
  documented with their mathematical convention and asymptotic complexity.
- Shared exact integer predicates and number-theoretic helpers, such as
  squarefreeness tests reused by endomorphism-side arithmetic, may also live
  here when they are domain-agnostic and more than one mathematical subtree
  can reasonably consume them.
- Small exact prime-power infrastructure also fits here when it is shared
  across domains: examples include exponentiation-by-squaring on `BigUint`,
  cached ladders `1, ℓ, ..., ℓ^e`, and normalized prime-power factorizations
  that validate only integer structure rather than curve semantics.
- The same is true for tiny exact integer-arithmetic helpers such as
  `gcd`/`lcm` on `BigUint`: if several domains may need them, prefer one small
  canonical implementation here instead of duplicating local versions.
- Avoid adding shared helpers whose main purpose is to downcast mathematical
  scalars, orders, or degrees into fixed-width integers. Prefer `BigUint` /
  `BigInt` for exact public numerical surfaces, and keep any unavoidable
  memory-sized index conversion local to the algorithm that materializes the
  corresponding table or loop.
- The same policy applies to small exact `usize` helpers such as Euclidean
  `gcd`/`lcm` and quotient families derived from distinct prime divisors when
  they support curve/order algorithms but do not depend on curve semantics.
- If callers repeatedly need both “factor this integer” and “normalize/validate
  the resulting prime powers”, prefer exposing that as one method on the
  normalized factorization type itself rather than leaving a second free helper
  in a consumer module.
- The same applies to tiny exact `p`-adic / `\ell`-adic integer helpers such
  as valuations `v_\ell(n)`, provided the API keeps the arithmetic surface
  small and validates any “prime input” precondition honestly.
- For Hensel-style lifting helpers, keep the first shared surface exact,
  crate-internal, and explicit about whether it only handles the simple-root
  case `f'(x) != 0 mod p`; prefer trace/report value objects that record the
  chosen correction digits before broadening to singular or domain-specific
  variants.
- For fast Hensel routes that double precision, record both the source and
  target levels of each step instead of pretending the route is a sequence of
  adjacent `k -> k + 1` lifts.
- For square roots modulo `2^e`, keep the two-adic route separate from the
  odd-prime simple-root route, since `2x` is never a unit modulo `2`; prefer
  explicit bit-lifting and returning all canonical roots.
- For square roots modulo a general integer `m`, prefer the canonical route
  “factor `m = Π pᵢ^eᵢ`, solve each prime-power component, recombine all local
  roots by CRT” under `numerics`, instead of duplicating that composition in
  curve-side algorithms.
- When testing square roots modulo general `m`, prefer the shared
  `proptest_support::numerics` generators with a brute-force oracle for small
  moduli before adding one-off local random-case builders.
- Keep Cornacchia-style Diophantine algorithms under their own `numerics`
  module rather than inside Hensel: they may consume modular square roots, but
  their core responsibility is solving equations such as `x² + d y² = m`.
- For Cornacchia over all roots modulo `m`, keep the candidate surface distinct
  from the primitive-solution surface: when `m` is not square-free, not every
  candidate solution need be primitive.
- For public representation questions such as `m = x² + d y²`, prefer the
  `quadratic_forms` value-object layer over exposing callers directly to the
  Cornacchia engine; Cornacchia remains the current primitive-representation
  implementation strategy.
- If an algorithm needs an exact integer square-root check, reuse the
  crate-internal helpers in `number_theory` rather than adding a local binary
  search or conflating the check with modular square roots.
- If an exact numerical algorithm is currently only a tested prototype with no
  production consumer, prefer gating that module or its one-off polynomial
  helpers behind `#[cfg(test)]` instead of letting normal builds accumulate
  dead-code warnings. Promote it back into crate-internal builds when a real
  caller appears.
- If a shared arithmetic helper enumerates a finite integer set such as the
  positive divisors of `n`, prefer documenting the chosen ordering convention
  explicitly and keeping the surface small and exact instead of wrapping it in
  heavier number-theory abstractions.
- Shared exact linear recurrences, such as order-2 recurrences evaluated
  either termwise or through one prefix pass, also fit here when several
  domains can reuse the same companion-matrix viewpoint.
- When such a recurrence needs both a fast isolated-term API and a batched
  prefix API, prefer one small value object that exposes both surfaces and
  documents clearly which one is `Θ(log n)` and which one is `Θ(N)`.
- When a classical arithmetic function needs more than one honest evaluation
  strategy, such as divisor-power sums `σ_k(n)`, it is acceptable to keep
  multiple algorithms side by side in `numerics`:
  a literal reference implementation, a better single-input algorithm, and a
  batched `1..=N` algorithm for coefficient-table generation.
- When a shared exact sequence has competing sign conventions, as Bernoulli
  numbers do at `B₁`, document the chosen convention explicitly and normalize
  the implementation to that convention instead of leaving callers to infer it
  from tests.
- For exact arithmetic helpers that feed symbolic coefficient formulas,
  prefer arbitrary-precision integer and rational backends (`BigInt`,
  `BigRational`) over fixed-width integers whenever growth can naturally exceed
  educational toy examples.
- If a domain algorithm needs a minimal integer denominator-clearing scale for
  an exact rational expression, prefer a shared `BigRational` helper here that
  documents the arithmetic convention, for example the least `u` such that
  `uᵏq ∈ ℤ`, instead of factoring denominators locally in the domain module.
  For simple rational-integrality checks, prefer `BigRational::is_integer()`
  directly over wrapping the denominator check in a project helper.
- If several domains need exact integer gcd/lcm behavior, keep the
  `BigInt`/`BigUint` helpers here and reexport them crate-privately instead of
  defining local Euclidean loops in polynomial or curve modules.
- If several exact algorithms need to cross from a mathematically positive
  `BigInt` into a `BigUint`, use the crate-internal
  `positive_bigint_to_biguint` helper instead of defining local conversion
  shims at each call site.
- Compatible Chinese-remainder reconstruction for non-coprime moduli belongs
  here when exact algebraic consumers need it, but keep it `pub(crate)` until
  there is a clear external API story beyond the existing public coprime CRT
  surface.
- When a shared exact helper depends on the external integer-factorization
  backend, still write its rustdoc cost in `Θ(...)` notation, but prefer a
  readable coarse term such as `factor(n)` over a crowded parameter list.
- Brute-force test oracles may use memory-sized loop bounds when the sampled
  domain is intentionally tiny, but compute residues and expected arithmetic
  values with `BigInt`/`BigUint` so the oracle does not normalize through
  fragile fixed-width casts.
- If one consumer needs all values up to a truncation bound, prefer one
  documented batched routine over repeated calls to a single-input helper, and
  say so explicitly in the consuming rustdocs.
- If a type is only meaningful inside one domain and has no shared numerical
  role, keep it local to that domain instead.
- Reexport from consumer modules when that improves ergonomics, but keep the
  canonical definition here when multiple domains depend on it.
- Keep the public root surface of `numerics` deliberately small:
  - tolerances, reusable complex-path geometry, validated quadrature objects,
    and a few standalone exact helpers with obvious end-user value are good
    public candidates
  - reference-only algorithms, backend glue, internal comparison wrappers,
    exact factorization scaffolding, and crate-internal arithmetic helpers
    should prefer `pub(crate)` even when several internal domains reuse them

## Testing expectations

- Test documented presets directly.
- Test constructor storage behavior directly.
- Keep expected constants explicit in the tests.
- When a tolerance helper implements a mixed absolute/relative comparison,
  test both a near-zero case and a scale-sensitive large-magnitude case.
- For reusable complex-path helpers, test the actual parameterization formulas
  directly, including the `s / (1 - s)` compactification used for rays toward
  infinity.
- For shared quadrature helpers, test both exact low-degree reference cases
  and any documented evaluation-order guarantees that downstream callers may
  rely on.

## Documentation expectations

- State clearly whether a constructor validates inputs or only packages them.
- Say whether a type is exact policy, approximate policy, or experimental
  scaffolding.
- Document the comparison rule explicitly when a helper uses more than a plain
  absolute epsilon.
- For exact-sequence algorithms such as Bernoulli-number generation, document
  the mathematical recurrence or tableau directly and state complexity in
  `Θ(...)` notation, making clear whether the estimate is in arithmetic
  operations, memory, or bit-complexity.
