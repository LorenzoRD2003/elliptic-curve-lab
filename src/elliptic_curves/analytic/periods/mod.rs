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
//!
//! Current internal architecture:
//!
//! - `root_recovery/` recovers the roots of the Weierstrass cubic
//! - `roots/` owns validated cubic-root triples plus their local
//!   classification and root-to-Legendre reductions
//! - `legendre/` chooses and analyzes one Legendre normalization
//! - `elliptic_integral/` evaluates `K(λ)` and related AGM-side data
//! - `period_basis/` transports those normalized periods back to the original
//!   curve and packages the resulting lattice reports
mod config;
pub mod elliptic_integral;
pub mod legendre;
mod metadata;
pub mod period_basis;
pub mod root_recovery;
pub mod roots;

pub use config::PeriodRecoveryConfig;
pub use legendre::{LegendreParameter, LegendreReduction, LegendreReductionReport};
pub use metadata::{NumericalRecoveryMetadata, PeriodRecoveryMethod, PeriodRecoveryStatus};
pub use period_basis::{
    CanonicalTauRecoveryReport, PeriodBasisRecoveryReport, RecoveredPeriodBasis,
    RecoveredPeriodBasisReport, TauRecoveryReport,
};
pub use root_recovery::CubicRootRecoveryReport;
pub use roots::WeierstrassCubicRoots;

pub(crate) use elliptic_integral::CompleteEllipticIntegralKApprox;
pub(crate) use period_basis::PeriodLatticeApprox;
pub(crate) use roots::{CubicRootConfiguration, CubicRootConfigurationReport, CubicRootSeparation};

#[cfg(test)]
mod tests;
