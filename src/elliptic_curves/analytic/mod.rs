//! Educational scaffolding for the future complex-analytic elliptic-curve
//! milestone.

pub mod comparison;
pub mod eisenstein;
pub mod elliptic_functions;
pub mod errors;
pub mod fundamental_domain;
pub mod invariants;
pub mod lab_report;
pub mod lattice;
pub mod modular_action;
pub mod periods;
pub mod q_expansion;
pub mod torsion;
pub mod torus_point;
pub mod truncation;
pub mod upper_half_plane;
pub mod weierstrass_model;
pub mod zeta;

pub use crate::numerics::ApproxTolerance;
pub use crate::numerics::tolerance;
pub use comparison::{
    ComplexApproxComparison, ComplexDifferenceReport, HasComplexApproxComparison,
};
pub use eisenstein::{
    EisensteinSumApprox, TruncationConvergenceReport, compare_eisenstein_truncations,
    eisenstein_sum, g4_sum, g6_sum,
};
pub use elliptic_functions::{
    EllipticFunctionApproximation, EllipticFunctionTruncation, HasPoleDistance, WeierstrassPApprox,
    WeierstrassPDerivativeApprox, evaluate_truncated_elliptic_function, weierstrass_p,
    weierstrass_p_derivative,
};
pub use errors::AnalyticCurveError;
pub use fundamental_domain::{
    FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
    FundamentalDomainReductionStep, FundamentalDomainReductionStepReason,
    is_in_standard_fundamental_domain, reduce_tau_to_standard_fundamental_domain,
};
pub use invariants::{
    AnalyticInvariants, analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
    analytic_invariants_from_tau, analytic_j_invariant,
};
pub use lab_report::{
    AnalyticShortWeierstrassModel, ComplexAnalyticCurveLabReport, HasAnalyticLatticeContext,
    SpecialJKind, SpecialTauKind, UniformizationExperimentReport,
};
pub use lattice::{
    ComplexLattice, ComplexTorusPoint, FundamentalParallelogramCoordinate, LatticeIndexPoint,
};
pub use modular_action::{ModularInvarianceReport, ModularMatrix, verify_j_modular_invariance};
pub use q_expansion::{
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
    JInvariantComparisonReport, JInvariantQExpansion, JInvariantQExpansionApprox,
    ModularQExpansionApproximation, ModularQExpansionCoefficients, ModularQExpansionFamily,
    ModularQParameter, QExpansionTruncation, compare_j_from_eisenstein_and_q_expansion,
};
pub use torsion::{
    AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
    AnalyticEvenDivisionPolynomialReport, AnalyticOddDivisionPolynomialReport,
    AnalyticTorsionPointApprox, EvenDivisionPolynomialVanishingBranch, TorusTorsionIndex,
    TorusTorsionPoint, compare_analytic_torsion_with_division_polynomial,
    compare_primitive_analytic_torsion_with_division_polynomial,
    map_primitive_torus_torsion_to_curve, map_torus_torsion_to_curve,
    primitive_torus_n_torsion_points, torus_n_torsion_points,
};
pub use torus_point::{
    TorusToCurveMapResult, TorusToCurveValues, WeierstrassDifferentialEquationReport,
    WeierstrassDifferentialEquationStatus, map_fundamental_point_to_curve,
    map_torus_point_to_curve, verify_weierstrass_differential_equation,
};
pub use truncation::LatticeSumTruncation;
pub use upper_half_plane::UpperHalfPlanePoint;
pub use weierstrass_model::{
    AnalyticCurveMembershipReport, AnalyticCurvePoint, AnalyticWeierstrassCurve,
};
pub use zeta::{WeierstrassZetaApprox, weierstrass_zeta};
