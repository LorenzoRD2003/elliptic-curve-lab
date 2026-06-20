//! General Weierstrass curves introduced in staged milestones.
//!
//! The long-term target is the affine model
//!
//! `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`.
//!
//! The current milestone already supports:
//!
//! - validated curve descriptors and classical invariants
//! - affine membership and the minimal `CurveModel` capability
//! - explicit conversions to and from short Weierstrass form
//!
//! while still deferring the general-model group law and the larger executable
//! stack that the mature short-Weierstrass family already owns.

mod display;
mod invariants;
mod membership;
mod model_traits;
mod reduction;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::GeneralWeierstrassCurve;
