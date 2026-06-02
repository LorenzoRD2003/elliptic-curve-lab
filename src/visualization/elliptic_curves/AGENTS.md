# AGENTS.md for `src/visualization/elliptic_curves`

## Module mission

This subtree explains curve models, points, and group-law behavior as readable
plain text.

Its output should help a learner see both:

- the algebraic object being displayed
- the specific geometric or arithmetic rule being used

## Scope guidance

- Prefer helpers that explain short-Weierstrass curves concretely before
  generalizing.
- Keep group-law explanations explicit about cases such as identity, inverse
  points, secants, and tangents.
- For small finite groups, it is appropriate to list points and describe point
  orders, but say directly when the method is exhaustive or based on repeated
  addition.
- If a helper depends on a backend-specific capability, such as square roots
  or enumeration, surface that honestly in the explanation text.
- For short-Weierstrass isomorphism and twist explanations, keep the
  distinction explicit between:
  same `j`-invariant = isomorphic over an algebraic closure,
  versus
  isomorphic over the base field = an actual scaling witness was found in the
  current field.
- For milestone-7 division-polynomial helpers, keep the distinction explicit
  between:
  - rational `x`-candidates
  - lifted rational points
  - points satisfying `[n]P = O`
  - points of exact order `n`
  - comparison against exhaustive enumeration

## Formatting guidance

- Prefer compact equation strings such as `y^2 = x^3 + ax + b`.
- Use `O` for the point at infinity in compact output.
- In richer explanations, show the important intermediate values, such as the
  left and right sides of the curve equation or the slope used in point
  addition.
- For milestone-7 explanations, include the polynomial shape explicitly:
  `polinomio en x` versus `y` times a polynomial in `x`.
- Avoid decorative ASCII art unless it materially helps understanding.

## Testing expectations

- Test the compact formatters separately from the richer descriptions.
- Test important phrases, mathematical quantities, and honest caveats.
- Include at least one finite-field example where point addition, point
  listing, and point order are all explained consistently.
- Include at least one division-polynomial example where the explanation
  mentions roots, lifted points, exact-order torsion, and the final
  enumeration comparison.
