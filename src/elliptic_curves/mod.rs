//! Elliptic-curve scaffolding.
//!
//! The crate currently exposes a mature short-Weierstrass stack together with
//! a staged general-Weierstrass model that will grow incrementally.

pub mod affine;
pub mod analytic;
pub mod endomorphisms;
pub mod error;
pub mod frobenius;
mod group_algorithms;
pub mod models;
pub mod projective;
pub use crate::elliptic_curves::models::general_weierstrass;
pub use crate::elliptic_curves::models::montgomery;
pub use crate::elliptic_curves::models::short_weierstrass;
pub use crate::elliptic_curves::models::twisted_edwards;

pub mod traits {
    pub use crate::elliptic_curves::models::traits::*;
}

pub use affine::AffinePoint;
pub use error::CurveError;
pub use models::general_weierstrass::GeneralWeierstrassCurve;
pub use models::montgomery::MontgomeryCurve;
pub use models::short_weierstrass::ShortWeierstrassCurve;
pub use models::twisted_edwards::{ExtendedTwistedEdwardsPoint, TwistedEdwardsCurve};
pub use projective::{CoordinateOperationCost, ProjectivePoint};
