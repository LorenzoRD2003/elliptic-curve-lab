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
//!
//! The subfiles reflect that pipeline:
//!
//! - `q_parameter.rs` validates `τ ↦ q`.
//! - `truncation.rs` validates finite truncation policies.
//! - `coefficients.rs` stores exact truncated coefficient tables.
//! - `eisenstein_series.rs` and `j_invariant.rs` provide the concrete
//!   expansion families.
//! - `comparison.rs` compares two approximation routes to the same modular
//!   quantity.

mod coefficients;
mod comparison;
mod eisenstein_series;
mod family;
mod j_invariant;
mod q_parameter;
mod truncation;

#[cfg(test)]
mod tests;

pub use comparison::JInvariantComparisonReport;
pub use eisenstein_series::{
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
};
pub use j_invariant::{JInvariantQExpansion, JInvariantQExpansionApprox};
pub use q_parameter::ModularQParameter;
pub use truncation::QExpansionTruncation;

pub(crate) use coefficients::ModularQExpansionCoefficients;
pub(crate) use family::{ModularQExpansionApproximation, ModularQExpansionFamily};
