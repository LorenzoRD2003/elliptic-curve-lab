//! Reduction of upper-half-plane parameters to the standard modular domain.
//!
//! This namespace packages the report types and owner methods for choosing a
//! more canonical representative of the modular class of `τ`.

mod api;
mod types;

#[cfg(test)]
mod tests;

pub use types::{
    FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
    FundamentalDomainReductionStep, FundamentalDomainReductionStepReason,
};

#[cfg(test)]
use crate::{
    elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint},
    numerics::ApproxTolerance,
};

#[cfg(test)]
pub(crate) fn reduce_tau_to_standard_fundamental_domain(
    tau: UpperHalfPlanePoint,
    max_steps: usize,
) -> Result<FundamentalDomainReductionReport, AnalyticCurveError> {
    tau.reduce_to_standard_fundamental_domain(max_steps)
}

#[cfg(test)]
pub(crate) fn is_in_standard_fundamental_domain(
    tau: &UpperHalfPlanePoint,
    tolerance: ApproxTolerance,
) -> bool {
    tau.is_in_standard_fundamental_domain(tolerance)
}
