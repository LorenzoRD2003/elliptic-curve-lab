//! Rational torsion for short-Weierstrass curves over `Q`.
//!
//! The current route is exact and intentionally classical. For a curve
//! `E/Q : y² = x³ + ax + b`, [`ShortWeierstrassCurve::rational_torsion`]
//! first transports the model to an integral short-Weierstrass companion when
//! denominators require it, then applies the Lutz-Nagell finite search. Each
//! candidate is verified by exact scalar multiplication against the finite list
//! of point orders allowed by Mazur's theorem.
//!
//! The returned [`RationalTorsionReport`] records the source curve, the integral
//! companion, the scaling factor, the Mazur-shape classification, every
//! certified rational torsion point, and the number of Lutz-Nagell candidates
//! checked. Later stages will add the reduction-mod-`p` and Hensel route from
//! the problem-set exercise as a second comparable strategy.
//!
//! Internal helpers for integral-model construction, candidate enumeration, and
//! verification remain crate-private; external callers should start from the
//! curve method and inspect the report.

mod classification;
mod enumeration;
mod error;
mod integral_model;
mod mazur;
mod report;
mod verification;

#[cfg(test)]
mod tests;

pub use classification::{RationalTorsionGroup, RationalTorsionGroupShape};
pub use error::RationalTorsionError;
pub use report::RationalTorsionReport;
