//! Educational inverse-uniformization helpers for analytic elliptic curves.
//!
//! This module collects the two current directions that start from the
//! torus-side recovery data and go back toward the curve:
//!
//! - validating a recovered `τ` or recovered lattice against curve-side
//!   invariants
//! - scaffolding the point-level Abel-Jacobi inverse map `(x, y) -> z mod Λ`

mod abel_jacobi;
mod validation;

#[cfg(test)]
mod tests;

pub use abel_jacobi::{
    AbelJacobiConfig, AbelJacobiContourReport, AbelJacobiInitialBranchChoice,
    AbelJacobiIntegralApprox, AbelJacobiIntegralDecomposition, AbelJacobiIntegralNumerics,
    AbelJacobiPointRecoveryReport, AbelJacobiRecoveryMetadata, AbelJacobiRecoveryStatus,
    AbelJacobiRoundtripValidationReport, AbelJacobiValidationPolicy,
    InverseUniformizationPointRecoveryReport, LegendreContourStrategy,
    PointRoundTripValidationConfig, PointRoundTripValidationReport,
    approximate_abel_jacobi_integral, recover_torus_point_from_curve_point,
    recover_torus_point_from_curve_point_with_periods,
    validate_point_inverse_uniformization_roundtrip,
    validate_point_inverse_uniformization_roundtrip_with_periods,
};
pub use validation::{
    InvariantRecoveryInterpretation, InvariantRecoveryValidationReport,
    InverseUniformizationJValidationReport, validate_canonical_tau_recovery_by_j_invariant,
    validate_recovered_lattice_invariants, validate_recovered_tau_by_j_invariant,
    validate_tau_recovery_report_by_j_invariant,
};
