//! Complex-analytic companions to the algebraic elliptic-curve modules.
//!
//! This subtree studies an elliptic curve over `ℂ` from the uniformization
//! viewpoint: one starts from a lattice `Λ ⊂ ℂ`, forms the torus `ℂ / Λ`,
//! and recovers the cubic curve through Weierstrass functions, modular
//! invariants, period recovery, and inverse-uniformization experiments.
//!
//! The namespace is intentionally staged.
//!
//! - `lattice/` owns validated lattice data and torus-side coordinates.
//! - `elliptic_functions/` evaluates truncated `℘` and `℘′` sums.
//! - `uniformization/` maps torus representatives forward to the cubic.
//! - `periods/` goes in the opposite direction: roots, Legendre reduction,
//!   complete elliptic integrals, and recovered periods.
//! - `inverse_uniformization/` validates recovered `τ`/lattice data and
//!   approximates Abel-Jacobi recovery back to torus classes.
//! - `modular_action/` and `q_expansion/` organize the modular side of the
//!   story.
//! - `torsion/` compares analytic torus torsion with algebraic torsion data.
//!
//! The current code is educational and numerically approximate. The goal is
//! not production numerical analysis, but a transparent bridge between the
//! classical analytic formulas and the crate's algebraic elliptic-curve APIs.
pub mod eisenstein;
pub mod elliptic_functions;
pub mod errors;
pub mod fundamental_domain;
pub mod invariants;
pub mod inverse_uniformization;
pub mod lattice;
pub mod modular_action;
pub mod periods;
pub mod q_expansion;
pub mod torsion;
pub mod uniformization;
pub mod upper_half_plane;
pub mod weierstrass_model;
pub mod zeta;

pub use crate::numerics::ApproxTolerance;
pub use elliptic_functions::EllipticFunctionTruncation;
pub use errors::AnalyticCurveError;
pub use lattice::{
    ComplexLattice, ComplexTorusPoint, FundamentalParallelogramCoordinate, LatticeIndexPoint,
    LatticeSumTruncation,
};
pub use upper_half_plane::UpperHalfPlanePoint;
pub use weierstrass_model::{AnalyticCurvePoint, AnalyticWeierstrassCurve};

#[allow(unused_imports)]
pub(crate) use eisenstein::{EisensteinSumApprox, TruncationConvergenceReport};
#[allow(unused_imports)]
pub(crate) use elliptic_functions::{WeierstrassPApprox, WeierstrassPDerivativeApprox};
#[allow(unused_imports)]
pub(crate) use fundamental_domain::{
    FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
    FundamentalDomainReductionStep, FundamentalDomainReductionStepReason,
};
pub(crate) use invariants::AnalyticInvariants;
#[allow(unused_imports)]
pub(crate) use inverse_uniformization::{
    AbelJacobiConfig, AbelJacobiPointRecoveryReport, AbelJacobiRecoveryMetadata,
    InvariantRecoveryInterpretation, InvariantRecoveryValidationReport,
    InverseUniformizationJValidationReport, PointRoundTripValidationConfig,
    PointRoundTripValidationReport,
};
pub(crate) use lattice::ComplexModuloLatticeComparison;
#[allow(unused_imports)]
pub(crate) use modular_action::{ModularInvarianceReport, ModularMatrix};
pub(crate) use periods::elliptic_integral::LegendrePeriodIntegralReport;
#[allow(unused_imports)]
pub(crate) use periods::legendre::{
    LegendreOrbitElementKind, LegendreParameterConditioning, LegendreParameterOrbit,
};
#[allow(unused_imports)]
pub(crate) use periods::period_basis::{CurvePeriodLatticeComparisonReport, PeriodLatticeApprox};
#[allow(unused_imports)]
pub(crate) use periods::{
    CanonicalTauRecoveryReport, CubicRootRecoveryReport, LegendreParameter, LegendreReduction,
    LegendreReductionReport, NumericalRecoveryMetadata, PeriodBasisRecoveryReport,
    PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryStatus, RecoveredPeriodBasis,
    RecoveredPeriodBasisReport, TauRecoveryReport, WeierstrassCubicRoots,
};
#[allow(unused_imports)]
pub(crate) use periods::{
    CubicRootConfiguration, CubicRootConfigurationReport, CubicRootSeparation,
};
#[allow(unused_imports)]
pub(crate) use q_expansion::{JInvariantComparisonReport, ModularQParameter, QExpansionTruncation};
#[allow(unused_imports)]
pub(crate) use torsion::{AnalyticDivisionPolynomialComparisonCase, AnalyticTorsionPointApprox};
#[allow(unused_imports)]
pub(crate) use uniformization::{
    TorusToCurveMapResult, TorusToCurveValues, WeierstrassDifferentialEquationReport,
    WeierstrassDifferentialEquationStatus,
};
pub(crate) use weierstrass_model::AnalyticCurveMembershipReport;

#[cfg(test)]
pub(crate) use fundamental_domain::{
    is_in_standard_fundamental_domain, reduce_tau_to_standard_fundamental_domain,
};
#[cfg(test)]
pub(crate) use inverse_uniformization::AbelJacobiValidationPolicy;
#[cfg(test)]
pub(crate) use periods::elliptic_integral::{
    CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata, ComplexAgmBranchChoice,
    ComplexAgmConfig, ComplexAgmStatus, complex_agm, complex_agm_trace,
};
#[cfg(test)]
#[allow(unused_imports)]
pub(crate) use q_expansion::JInvariantQExpansion;
