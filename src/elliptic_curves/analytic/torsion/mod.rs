//! Analytic torsion on `ℂ / Λ` and its comparison with algebraic torsion.
//!
//! The analytic identification
//! `E[n] ≅ (1/n)Λ / Λ`
//! provides a clean bridge from torus-side torsion indices to the algebraic
//! `n`-torsion on the Weierstrass cubic.
//! This module packages that bridge for the small educational experiments in
//! the crate.
//!
//! The current split is:
//!
//! - `types.rs` for validated torus torsion indices and high-level comparison
//!   cases.
//! - `torus.rs` for torus-side torsion sampling and exact-order checks.
//! - `curve_map.rs` for the map from torus torsion to curve points.
//! - `division_polynomial.rs` for comparison against the algebraic
//!   division-polynomial side.
//!
//! The public API emphasizes explicit cases and reports so the user can see
//! when a torsion class maps to infinity, when an even-index subtlety forces
//! `y ≈ 0`, and how the analytic and algebraic viewpoints line up.
mod curve_map;
mod division_polynomial;
#[cfg(test)]
mod tests;
mod torus;
mod types;

pub use types::{
    AnalyticDivisionPolynomialComparisonCase, AnalyticTorsionPointApprox, TorusTorsionIndex,
    TorusTorsionPoint,
};

pub(crate) use types::{
    AnalyticDivisionPolynomialComparisonStatus, AnalyticEvenDivisionPolynomialReport,
    AnalyticOddDivisionPolynomialReport, EvenDivisionPolynomialVanishingBranch,
};
