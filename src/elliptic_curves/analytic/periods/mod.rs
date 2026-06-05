//! Approximate period-lattice metadata for future analytic period recovery.
//!
//! For a non-singular analytic Weierstrass model
//! `y² = 4x³ - g₂x - g₃`, one expects a complex lattice
//! `Λ = ℤω₁ + ℤω₂` whose associated torus `ℂ / Λ` uniformizes the curve.
//! The ordered basis is not unique: replacing `(ω₁, ω₂)` by another
//! positively oriented `SL₂(ℤ)`-equivalent basis describes the same lattice.
//!
//! This module therefore starts with small metadata objects rather than a
//! premature recovery algorithm. The current surface is meant to support:
//!
//! - storing one chosen approximate period basis
//! - recording the corresponding modulus `τ = ω₂ / ω₁`
//! - comparing the `j`-invariant implied by that recovered lattice against
//!   the original curve-side `j`
mod agm;
mod classification;
mod config;
mod elliptic_integral;
mod lattice;
mod legendre;
mod metadata;
mod recovery;
mod report;
mod roots;

pub use agm::{
    ComplexAgmBranchChoice, ComplexAgmConfig, ComplexAgmIteration, ComplexAgmResult,
    ComplexAgmStatus, ComplexAgmTrace, complex_agm, complex_agm_trace,
};
pub use classification::{
    CubicRootConfiguration, CubicRootConfigurationReport, CubicRootSeparation,
    classify_cubic_root_configuration, cubic_root_configuration_report,
};
pub use config::PeriodRecoveryConfig;
pub use elliptic_integral::{
    CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata,
    LegendrePeriodIntegralReport, complementary_complete_elliptic_integral_k_from_lambda,
    complementary_complete_elliptic_integral_k_from_m, complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_m, legendre_period_integral_report,
};
pub use lattice::PeriodLatticeApprox;
pub use legendre::{
    LegendreOrbitElement, LegendreOrbitElementKind, LegendreParameter,
    LegendreParameterConditioning, LegendreParameterOrbit, LegendreReduction,
    LegendreReductionReport, classify_legendre_parameter_conditioning, legendre_reduction_report,
};
pub use metadata::{NumericalRecoveryMetadata, PeriodRecoveryMethod, PeriodRecoveryStatus};
pub use recovery::{
    CubicRootRecoveryReport, recover_weierstrass_cubic_roots,
    recover_weierstrass_cubic_roots_from_invariants, recover_weierstrass_cubic_roots_with_report,
};
pub use report::PeriodRecoveryReport;
pub use roots::WeierstrassCubicRoots;

#[cfg(test)]
mod tests;
