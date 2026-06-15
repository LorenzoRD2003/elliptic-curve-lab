//! Short-Weierstrass division polynomials and rational torsion helpers.
//!
//! - construct the honest short-Weierstrass shape `ψ_n ∈ F[x]` or `ψ_n ∈ yF[x]`
//! - evaluate it at affine points
//! - recover small rational torsion through finite-field enumeration
//! - compare that recovery against exhaustive group enumeration

mod comparison;
mod construction;
mod criterion_dispatch;
mod error;
mod evaluation;
mod torsion_report;
mod torsion_search;
mod types;

#[cfg(test)]
mod tests;

pub use error::DivisionPolynomialError;
pub use torsion_report::TorsionComparisonReport;
pub use types::DivisionPolynomialForm;
