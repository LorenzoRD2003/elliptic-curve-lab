//! Twisted Edwards curves introduced in staged milestones.
//!
//! The current milestone supports:
//!
//! - validated curve descriptors for `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2`
//! - classical invariants
//! - affine membership, `CurveModel`, `AffineCurveModel`, `HasJInvariant`,
//!   and `LiftXCoordinate`
//! - total whole-curve conversion to and from the Montgomery family
//! - native affine group-law support
//! - compatibility with the shared small finite-group surfaces through
//!   `EnumerableCurveModel` and `FiniteGroupCurveModel`
//!
//! over fields of characteristic different from `2`.

mod group_law;
mod invariants;
mod membership;
mod model_traits;
mod reduction;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::TwistedEdwardsCurve;
pub use reduction::TwistedEdwardsBirationalMapError;
