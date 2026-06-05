//! Small building blocks for modular `q`-parameters and first `q`-expansions.
//!
//! The classical change of variable `q = e^{2π i τ}`
//! turns a point `τ` in the upper half-plane into a small complex parameter
//! inside the open unit disc, because `|q| = e^{-2π Im(τ)} < 1`
//! whenever `Im(τ) > 0`.
//!
//! This is the standard coordinate used in Fourier and `q`-expansions of
//! modular forms and modular functions near the cusp `i∞`.
//!
//! In particular:
//!
//! - `E₄` and `E₆` are holomorphic modular forms of weights `4` and `6`
//! - `j` is a modular function of weight `0`, not a holomorphic modular form
//!
//! So the shared abstraction in this module is intentionally phrased in terms
//! of neutral modular `q`-expansion families rather than only modular forms.
//!
//! In this crate we keep both the validated input `τ` and the derived value
//! `q` together, so later analytic routines can report both the geometric
//! modular parameter and the actual small complex quantity used in a numerical
//! expansion.

mod coefficients;
mod comparison;
mod eisenstein_series;
mod family;
mod j_invariant;
mod q_parameter;
mod truncation;

pub use coefficients::ModularQExpansionCoefficients;
pub use comparison::{JInvariantComparisonReport, compare_j_from_eisenstein_and_q_expansion};
pub use eisenstein_series::{
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
};
pub use family::{ModularQExpansionApproximation, ModularQExpansionFamily};
pub use j_invariant::{JInvariantQExpansion, JInvariantQExpansionApprox};
pub use q_parameter::ModularQParameter;
pub use truncation::QExpansionTruncation;

#[cfg(test)]
mod tests;
