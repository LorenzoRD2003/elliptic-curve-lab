//! Educational inverse-uniformization helpers for analytic elliptic curves.
//!
//! This module collects the two current directions that start from the
//! torus-side recovery data and go back toward the curve:
//!
//! - validating a recovered `τ` or recovered lattice against curve-side
//!   invariants
//! - scaffolding the point-level Abel-Jacobi inverse map `(x, y) -> z mod Λ`
//!
//! The current split is deliberately small:
//!
//! - `j_validation/` compares recovered modular data through the
//!   `j`-invariant.
//! - `lattice_invariants/` compares lattice-side invariants against the
//!   source analytic curve.
//! - `abel_jacobi/` owns the point-level inverse map and its validation
//!   reports.
//! - `validation_shared.rs` holds small glue reused by the validation routes.

pub mod abel_jacobi;
pub mod j_validation;
pub mod lattice_invariants;
mod validation_shared;

#[cfg(test)]
mod tests;

pub use abel_jacobi::{
    AbelJacobiBudgets, AbelJacobiConfig, AbelJacobiPointRecoveryReport, AbelJacobiRecoveryMetadata,
    AbelJacobiRecoveryStatus, AbelJacobiValidationPolicy, InverseUniformizationPointRecoveryReport,
    PointRoundTripValidationConfig, PointRoundTripValidationReport,
};
pub use j_validation::InverseUniformizationJValidationReport;
pub use lattice_invariants::{InvariantRecoveryInterpretation, InvariantRecoveryValidationReport};
