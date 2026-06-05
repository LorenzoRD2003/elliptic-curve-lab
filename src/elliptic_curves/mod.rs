//! Elliptic-curve scaffolding.

pub mod affine;
pub mod analytic;
pub mod division_polynomials;
pub mod error;
pub mod invariants;
pub mod isomorphisms;
pub mod short_weierstrass;
pub mod torsion;
pub mod traits;

pub use affine::AffinePoint;
pub use analytic::{
    AnalyticCurveError, AnalyticCurveMembershipReport, AnalyticCurvePoint,
    AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
    AnalyticEvenDivisionPolynomialReport, AnalyticInvariants, AnalyticOddDivisionPolynomialReport,
    AnalyticShortWeierstrassModel, AnalyticTorsionPointApprox, AnalyticWeierstrassCurve,
    ApproxTolerance, CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata,
    ComplexAgmBranchChoice, ComplexAgmConfig, ComplexAgmIteration, ComplexAgmResult,
    ComplexAgmStatus, ComplexAgmTrace, ComplexAnalyticCurveLabReport, ComplexApproxComparison,
    ComplexDifferenceReport, ComplexLattice, ComplexTorusPoint, CubicRootConfiguration,
    CubicRootConfigurationReport, CubicRootRecoveryReport, CubicRootSeparation,
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
    EisensteinSumApprox, EllipticFunctionApproximation, EllipticFunctionTruncation,
    EvenDivisionPolynomialVanishingBranch, FundamentalDomainReductionReport,
    FundamentalDomainReductionStatus, FundamentalDomainReductionStep,
    FundamentalDomainReductionStepReason, FundamentalParallelogramCoordinate,
    HasAnalyticLatticeContext, HasComplexApproxComparison, HasPoleDistance,
    JInvariantComparisonReport, JInvariantQExpansion, JInvariantQExpansionApprox,
    LatticeIndexPoint, LatticeSumTruncation, LegendreOrbitElement, LegendreOrbitElementKind,
    LegendreParameter, LegendreParameterConditioning, LegendreParameterOrbit,
    LegendrePeriodIntegralReport, LegendreReduction, LegendreReductionReport,
    ModularInvarianceReport, ModularMatrix, ModularQExpansionApproximation,
    ModularQExpansionCoefficients, ModularQExpansionFamily, ModularQParameter,
    NumericalRecoveryMetadata, PeriodLatticeApprox, PeriodRecoveryConfig, PeriodRecoveryMethod,
    PeriodRecoveryReport, PeriodRecoveryStatus, QExpansionTruncation, SpecialJKind, SpecialTauKind,
    TorusToCurveMapResult, TorusToCurveValues, TorusTorsionIndex, TorusTorsionPoint,
    TruncationConvergenceReport, UniformizationExperimentReport, UpperHalfPlanePoint,
    WeierstrassCubicRoots, WeierstrassDifferentialEquationReport,
    WeierstrassDifferentialEquationStatus, WeierstrassPApprox, WeierstrassPDerivativeApprox,
    WeierstrassZetaApprox, analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
    analytic_invariants_from_tau, analytic_j_invariant, classify_cubic_root_configuration,
    classify_legendre_parameter_conditioning, compare_analytic_torsion_with_division_polynomial,
    compare_eisenstein_truncations, compare_j_from_eisenstein_and_q_expansion,
    compare_primitive_analytic_torsion_with_division_polynomial, comparison as analytic_comparison,
    complementary_complete_elliptic_integral_k_from_lambda,
    complementary_complete_elliptic_integral_k_from_m, complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_m, complex_agm, complex_agm_trace,
    cubic_root_configuration_report, eisenstein as analytic_eisenstein, eisenstein_sum,
    elliptic_functions as analytic_elliptic_functions, errors as analytic_errors,
    evaluate_truncated_elliptic_function, fundamental_domain as analytic_fundamental_domain,
    g4_sum, g6_sum, invariants as analytic_invariants_module, is_in_standard_fundamental_domain,
    lattice as analytic_lattice, legendre_period_integral_report, legendre_reduction_report,
    map_fundamental_point_to_curve, map_primitive_torus_torsion_to_curve, map_torus_point_to_curve,
    map_torus_torsion_to_curve, modular_action as analytic_modular_action,
    periods as analytic_periods, primitive_torus_n_torsion_points,
    q_expansion as analytic_q_expansion, recover_weierstrass_cubic_roots,
    recover_weierstrass_cubic_roots_from_invariants, recover_weierstrass_cubic_roots_with_report,
    reduce_tau_to_standard_fundamental_domain, tolerance as analytic_tolerance,
    torsion as analytic_torsion, torus_n_torsion_points, torus_point as analytic_torus_point,
    truncation as analytic_truncation, upper_half_plane as analytic_upper_half_plane,
    verify_j_modular_invariance, verify_weierstrass_differential_equation,
    weierstrass_model as analytic_weierstrass_model, weierstrass_p, weierstrass_p_derivative,
    weierstrass_zeta,
};
pub use error::CurveError;
pub use invariants::HasJInvariant;
pub use isomorphisms::{
    CurveIsomorphism, CurveIsomorphismError, ShortWeierstrassIsomorphism,
    ShortWeierstrassQuadraticTwist, ShortWeierstrassTwist, TwistKind,
};
pub use short_weierstrass::ShortWeierstrassCurve;
pub use torsion::{point_has_exact_order, points_of_exact_order};
pub use traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure,
    FiniteGroupCurveModel, GroupCurveModel, LiftXCoordinate, PointIndexSampler,
};
