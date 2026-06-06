//! Educational scaffolding for the complex-analytic elliptic-curve layer.

pub mod comparison;
pub mod eisenstein;
pub mod elliptic_functions;
pub mod errors;
pub mod fundamental_domain;
pub mod invariants;
pub mod inverse_uniformization;
pub mod lab_report;
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
pub use inverse_uniformization::{
    AbelJacobiConfig, AbelJacobiContourReport, AbelJacobiInitialBranchChoice,
    AbelJacobiIntegralApprox, AbelJacobiIntegralDecomposition, AbelJacobiIntegralNumerics,
    AbelJacobiPointRecoveryReport, AbelJacobiRecoveryMetadata, AbelJacobiRecoveryStatus,
    AbelJacobiRoundtripValidationReport, AbelJacobiValidationPolicy,
    InvariantRecoveryInterpretation, InvariantRecoveryValidationReport,
    InverseUniformizationJValidationReport, InverseUniformizationPointRecoveryReport,
    LegendreContourStrategy, PointRoundTripValidationConfig, PointRoundTripValidationReport,
    approximate_abel_jacobi_integral, recover_torus_point_from_curve_point,
    recover_torus_point_from_curve_point_with_periods,
    validate_canonical_tau_recovery_by_j_invariant,
    validate_point_inverse_uniformization_roundtrip,
    validate_point_inverse_uniformization_roundtrip_with_periods,
    validate_recovered_lattice_invariants, validate_recovered_tau_by_j_invariant,
    validate_tau_recovery_report_by_j_invariant,
};
pub use lab_report::{
    AnalyticShortWeierstrassModel, ComplexAnalyticCurveLabReport, HasAnalyticLatticeContext,
    SpecialJKind, SpecialTauKind, UniformizationExperimentReport,
};
pub use lattice::{
    ComplexLattice, ComplexModuloLatticeComparison, ComplexTorusPoint,
    FundamentalParallelogramCoordinate, LatticeIndexPoint, LatticeSumTruncation,
};
pub use modular_action::{ModularInvarianceReport, ModularMatrix, verify_j_modular_invariance};
pub use periods::{
    CanonicalTauRecoveryReport, CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata,
    ComplexAgmBranchChoice, ComplexAgmConfig, ComplexAgmIteration, ComplexAgmResult,
    ComplexAgmStatus, ComplexAgmTrace, CubicRootConfiguration, CubicRootConfigurationReport,
    CubicRootRecoveryReport, CubicRootSeparation, LegendreOrbitElement, LegendreOrbitElementKind,
    LegendreParameter, LegendreParameterConditioning, LegendreParameterOrbit,
    LegendrePeriodIntegralReport, LegendreReduction, LegendreReductionReport,
    NumericalRecoveryMetadata, PeriodBasisRecoveryReport, PeriodLatticeApprox,
    PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryReport, PeriodRecoveryStatus,
    RecoveredPeriodBasis, RecoveredPeriodBasisReport, TauRecoveryReport, WeierstrassCubicRoots,
    classify_cubic_root_configuration, classify_legendre_parameter_conditioning,
    complementary_complete_elliptic_integral_k_from_lambda,
    complementary_complete_elliptic_integral_k_from_m, complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_m, complex_agm, complex_agm_trace,
    cubic_root_configuration_report, legendre_period_integral_report, legendre_reduction_report,
    recover_canonical_tau_from_curve, recover_period_basis,
    recover_period_basis_from_legendre_reduction, recover_tau_from_curve,
    recover_weierstrass_cubic_roots, recover_weierstrass_cubic_roots_from_invariants,
    recover_weierstrass_cubic_roots_with_report,
};
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
pub use uniformization::{
    TorusToCurveMapResult, TorusToCurveValues, WeierstrassDifferentialEquationReport,
    WeierstrassDifferentialEquationStatus, map_fundamental_point_to_curve,
    map_torus_point_to_curve, verify_weierstrass_differential_equation,
};
pub use upper_half_plane::UpperHalfPlanePoint;
pub use weierstrass_model::{
    AnalyticCurveMembershipReport, AnalyticCurvePoint, AnalyticWeierstrassCurve,
};
pub use zeta::{WeierstrassZetaApprox, weierstrass_zeta};
