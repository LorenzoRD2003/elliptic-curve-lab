//! General Weierstrass curves introduced in staged milestones.
//!
//! The long-term target is the affine model
//!
//! `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`.
//!
//! The current milestone already supports:
//!
//! - validated curve descriptors and classical invariants
//! - affine membership together with `CurveModel`, `AffineCurveModel`, and
//!   `HasJInvariant`
//! - finite-field `LiftXCoordinate` support across both odd characteristic and
//!   characteristic `2`
//! - `LiftXCoordinate` over characteristic-`0` backends such as `Q` and
//!   `ComplexApprox` through the same completing-square story
//! - explicit conversions to and from short Weierstrass form
//! - `GroupCurveModel` support through native affine formulas for the general
//!   model
//! - compatibility with `EnumerableCurveModel`, `FiniteGroupCurveModel`, and
//!   `FrobeniusTraceCurveModel` over small finite fields
//! - compatibility with explicit `IsogenyKernel` construction from general
//!   Weierstrass points in those same small finite settings
//!
//! while still deferring the projective-coordinate general-model group law and
//! the larger executable stack that the mature short-Weierstrass family already
//! owns.
//!
//! The `LiftXCoordinate` story is now fiber-oriented:
//!
//! - in odd characteristic, the `y`-quadratic is solved by completing the
//!   square
//! - in characteristic `2`, the `u = 0` branch uses the unique Frobenius
//!   square root and the `u != 0` branch reduces to an Artin-Schreier solve
//!
//! Rational extension fields are intentionally still outside that surface for
//! now: the generic extension-field backend over `Q` does not yet expose a
//! square-root capability, so the odd-characteristic completing-square route
//! cannot be closed honestly there.

mod display;
mod group_law;
mod invariants;
mod membership;
mod model_traits;
mod reduction;
mod type_definition;
mod y_fiber;

#[cfg(test)]
mod tests;

pub use type_definition::GeneralWeierstrassCurve;
