//! Rational torsion for short-Weierstrass curves over `Q`.
//!
//! For a curve `E/ℚ : y² = x³ + ax + b`, the public entry point is
//! [`ShortWeierstrassCurve::rational_torsion_by`]. Callers choose an explicit
//! [`RationalTorsionStrategy`] so the report records whether the computation
//! used Lutz-Nagell candidate enumeration or the good-reduction/Hensel route.
//!
//! The returned [`RationalTorsionReport`] records the source curve, the integral
//! companion, the scaling factor, the Mazur-shape classification, every
//! certified rational torsion point, and strategy-specific metadata.
//!
//! Roadmap: a later educational variant may replace the Mazur order bound by a
//! reduction-gcd bound. The intended route is to choose several good primes
//! `pᵢ`, compute `gcd(#E(𝔽_{pᵢ}))`, and verify divisors of that bound exactly.
//! That route is useful for explaining why rational torsion injects into good
//! reductions, but it is not currently the default because Mazur's theorem gives
//! a smaller and clearer bound for curves over `ℚ`.
//!
//! Internal helpers for integral-model construction, candidate enumeration, and
//! verification remain crate-private; external callers should start from
//! `rational_torsion_by(...)` and inspect the report.

mod classification;
mod enumeration;
mod error;
mod integral_model;
mod mazur;
mod reduction_mod_p;
mod report;
mod strategy;
mod verification;

#[cfg(test)]
mod tests;

pub use classification::{RationalTorsionGroup, RationalTorsionGroupShape};
pub use error::RationalTorsionError;
pub use report::RationalTorsionReport;
pub use strategy::RationalTorsionStrategy;
