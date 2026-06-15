//! Recovery and packaging of period bases for the analytic curve.
//!
//! Once the cubic has been reduced to Legendre form and the corresponding
//! complete elliptic integrals have been approximated, the remaining task is
//! to transport those normalized periods back to the original analytic
//! Weierstrass model and package the result in user-facing reports.
//!
//! The current split is:
//!
//! - `value.rs` for recovered period bases as validated lattices.
//! - `lattice.rs` for lightweight packaged period-lattice approximations.
//! - `basis_report.rs` for the main period-basis recovery reports.
//! - `tau_report.rs` and `canonical_tau_report.rs` for the `τ`-focused views.
//! - `validation.rs` for comparing a recovered lattice against the curve-side
//!   `j`-invariant.
//! - `recovery.rs` for the inherent recovery methods on
//!   `LegendreReduction` and `AnalyticWeierstrassCurve`.
//! - `internal.rs` for crate-private assembly helpers.
//!
//! This submodule is the point where Legendre-side normalized data becomes the
//! full lattice data used by the rest of the analytic APIs.
mod basis_report;
mod canonical_tau_report;
mod internal;
mod lattice;
mod recovery;
mod tau_report;
mod validation;
mod value;

pub use basis_report::{PeriodBasisRecoveryReport, RecoveredPeriodBasisReport};
pub use canonical_tau_report::CanonicalTauRecoveryReport;
pub use lattice::PeriodLatticeApprox;
pub use tau_report::TauRecoveryReport;
pub use validation::CurvePeriodLatticeComparisonReport;
pub use value::RecoveredPeriodBasis;
