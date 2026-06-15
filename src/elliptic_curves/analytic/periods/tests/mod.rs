//! Internal tests for the analytic period-recovery pipeline.
//!
//! The files here are organized by the same mathematical stages used by the
//! implementation:
//!
//! - `config.rs` for validated recovery policies,
//! - `roots_and_legendre.rs` for cubic roots and Legendre normalization,
//! - `elliptic_integral.rs` and `agm.rs` for the normalized period side,
//! - `period_basis.rs` and `metadata_and_lattice.rs` for transported periods
//!   and lattice-level reports,
//! - `root_recovery.rs` and `end_to_end.rs` for the full recovery pipeline.
//!
//! This module is test-only, but the split mirrors the production ownership
//! boundaries so regressions stay close to the relevant mathematical story.
use num_complex::Complex64;
use proptest::prelude::*;

pub(super) use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, CanonicalTauRecoveryReport,
    CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata, ComplexAgmBranchChoice,
    ComplexAgmConfig, ComplexAgmStatus, ComplexLattice, CubicRootConfiguration,
    CubicRootConfigurationReport, CubicRootRecoveryReport, CubicRootSeparation,
    LatticeSumTruncation, LegendreOrbitElementKind, LegendreParameter,
    LegendreParameterConditioning, LegendreReduction, LegendreReductionReport,
    NumericalRecoveryMetadata, PeriodBasisRecoveryReport, PeriodRecoveryConfig,
    PeriodRecoveryMethod, PeriodRecoveryStatus, RecoveredPeriodBasis, RecoveredPeriodBasisReport,
    TauRecoveryReport, UpperHalfPlanePoint, WeierstrassCubicRoots, complex_agm, complex_agm_trace,
    is_in_standard_fundamental_domain,
};
pub(super) use crate::fields::complex_approx::ComplexApprox;
pub(super) use crate::proptest_support::elliptic_curves::arb_stable_real_split_analytic_curve;

mod agm;
mod config;
mod elliptic_integral;
mod end_to_end;
mod metadata_and_lattice;
mod period_basis;
mod root_recovery;
mod roots_and_legendre;
