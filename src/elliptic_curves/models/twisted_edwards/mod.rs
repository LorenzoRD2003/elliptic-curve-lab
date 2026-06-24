//! Twisted Edwards curves introduced in staged milestones.
//!
//! The current milestone supports:
//!
//! - validated curve descriptors for
//!   `a x^2 + y^2 = 1 + d x^2 y^2`
//! - classical invariants
//! - affine membership, `CurveModel`, `AffineCurveModel`, `HasJInvariant`,
//!   and `LiftXCoordinate`
//!
//! over fields of characteristic different from `2`.

mod invariants;
mod membership;
mod model_traits;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::TwistedEdwardsCurve;
