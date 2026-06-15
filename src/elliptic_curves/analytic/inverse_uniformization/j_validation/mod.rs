//! Validation of recovered modular parameters through the `j`-invariant.
//!
//! This submodule owns the `j`-level validation story for inverse
//! uniformization:
//!
//! - `report.rs` stores the recovered lattice data and the resulting
//!   `j(Λ_τ)` versus `j(E)` comparison.
//! - `api.rs` attaches the main validation entry points to
//!   `AnalyticWeierstrassCurve`, `TauRecoveryReport`, and
//!   `CanonicalTauRecoveryReport`.

mod api;
mod report;

pub use report::InverseUniformizationJValidationReport;
