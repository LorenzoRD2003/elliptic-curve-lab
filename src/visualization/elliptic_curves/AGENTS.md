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

## Formatting guidance

- Prefer compact equation strings such as `y^2 = x^3 + ax + b`.
- Use `O` for the point at infinity in compact output.
- In richer explanations, show the important intermediate values, such as the
  left and right sides of the curve equation or the slope used in point
  addition.
- Avoid decorative ASCII art unless it materially helps understanding.

## Testing expectations

- Test the compact formatters separately from the richer descriptions.
- Test important phrases, mathematical quantities, and honest caveats.
- Include at least one finite-field example where point addition, point
  listing, and point order are all explained consistently.
