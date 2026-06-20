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
//! - explicit conversions to and from short Weierstrass form
//!
//! while still deferring the general-model group law and the larger executable
//! stack that the mature short-Weierstrass family already owns.
//!
//! In particular, this stage intentionally does not implement
//! `LiftXCoordinate`: the current trait models equations of the form
//! `y^2 = rhs(x)`, while the general Weierstrass equation determines `y`
//! through a quadratic equation whose honest generic solution would need extra
//! structure, especially in characteristic `2`.

mod display;
mod invariants;
mod membership;
mod model_traits;
mod reduction;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::GeneralWeierstrassCurve;
