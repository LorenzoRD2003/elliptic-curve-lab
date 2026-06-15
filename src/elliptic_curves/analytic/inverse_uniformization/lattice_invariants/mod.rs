//! Validation of recovered lattices through full analytic invariants.
//!
//! This submodule owns the scale-sensitive validation story:
//!
//! - `report.rs` stores the recovered invariant comparisons and the
//!   interpretation of whether the mismatch is only a homothety issue or a
//!   genuine modular inconsistency.
//! - `api.rs` attaches the main validation entry point to
//!   `AnalyticWeierstrassCurve`.

mod api;
mod report;

pub use report::{InvariantRecoveryInterpretation, InvariantRecoveryValidationReport};
