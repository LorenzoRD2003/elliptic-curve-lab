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
- If callers repeatedly need both “factor this integer” and “normalize/validate
  the resulting prime powers”, prefer exposing that as one method on the
  normalized factorization type itself rather than leaving a second free helper
  in a consumer module.
- The same applies to tiny exact `p`-adic / `\ell`-adic integer helpers such
  as valuations `v_\ell(n)`, provided the API keeps the arithmetic surface
  small and validates any “prime input” precondition honestly.
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
- If one consumer needs all values up to a truncation bound, prefer one
  documented batched routine over repeated calls to a single-input helper, and
  say so explicitly in the consuming rustdocs.
- If a type is only meaningful inside one domain and has no shared numerical
  role, keep it local to that domain instead.
- Reexport from consumer modules when that improves ergonomics, but keep the
  canonical definition here when multiple domains depend on it.

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
