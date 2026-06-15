//! Elliptic-curve scaffolding.

pub mod affine;
pub mod analytic;
pub mod endomorphisms;
pub mod error;
pub mod frobenius;
mod group_algorithms;
pub mod models;
pub use crate::elliptic_curves::models::short_weierstrass;

pub mod traits {
    pub use crate::elliptic_curves::models::traits::*;
}

pub use affine::AffinePoint;
pub use error::CurveError;
pub use models::short_weierstrass::ShortWeierstrassCurve;
