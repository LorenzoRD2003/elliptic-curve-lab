//! Elliptic-curve scaffolding.

pub mod affine;
pub mod analytic;
pub mod division_polynomials;
pub mod endomorphisms;
pub mod error;
pub mod frobenius;
pub mod function_fields;
pub mod invariants;
pub mod isomorphisms;
mod order_from_multiple;
pub mod short_weierstrass;
pub mod torsion;
pub mod traits;

pub use affine::AffinePoint;
pub use analytic::{
    AbelJacobiConfig, AnalyticCurveError, AnalyticCurvePoint, AnalyticInvariants,
    AnalyticShortWeierstrassModel, AnalyticWeierstrassCurve, ApproxTolerance,
    CompleteEllipticIntegralKApprox, ComplexLattice, ComplexModuloLatticeComparison,
    ComplexTorusPoint, CubicRootConfiguration, EisensteinSeriesQExpansion, EisensteinSeriesWeight,
    EllipticFunctionTruncation, FundamentalParallelogramCoordinate, LatticeIndexPoint,
    LatticeSumTruncation, LegendreParameter, LegendreReduction, ModularMatrix, ModularQParameter,
    PeriodRecoveryConfig, QExpansionTruncation, RecoveredPeriodBasis, TorusTorsionIndex,
    TorusTorsionPoint, UpperHalfPlanePoint, WeierstrassCubicRoots, analytic_discriminant,
    analytic_g2, analytic_g3, analytic_invariants, analytic_invariants_from_tau,
    analytic_j_invariant, approximate_abel_jacobi_integral,
    compare_analytic_torsion_with_division_polynomial,
    compare_primitive_analytic_torsion_with_division_polynomial,
    complementary_complete_elliptic_integral_k_from_lambda,
    complementary_complete_elliptic_integral_k_from_m, complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_m, eisenstein_sum, g4_sum, g6_sum,
    map_fundamental_point_to_curve, map_primitive_torus_torsion_to_curve, map_torus_point_to_curve,
    map_torus_torsion_to_curve, primitive_torus_n_torsion_points, recover_canonical_tau_from_curve,
    recover_period_basis, recover_period_basis_from_legendre_reduction, recover_tau_from_curve,
    recover_torus_point_from_curve_point, recover_torus_point_from_curve_point_with_periods,
    recover_weierstrass_cubic_roots, recover_weierstrass_cubic_roots_from_invariants,
    reduce_tau_to_standard_fundamental_domain, torus_n_torsion_points, weierstrass_p,
    weierstrass_p_derivative, weierstrass_zeta,
};
pub use endomorphisms::{
    DiscriminantSign, EndomorphismRingCandidateSet, EndomorphismRingLocalView,
    EndomorphismRingReport, ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError,
    QuadraticDiscriminant, QuadraticDiscriminantFactorization,
    QuadraticDiscriminantFactorizationError, QuadraticDiscriminantMod4,
    QuadraticOrderCoverRelation, QuadraticOrderIndexError, VolcanoEndomorphismLevelCandidate,
};
pub use error::CurveError;
pub use frobenius::{
    AbsoluteFrobenius, CharacterSumPointCount, FrobeniusCharacteristicEquationCheck,
    FrobeniusCharacteristicEquationExhaustiveReport, FrobeniusCharacteristicPolynomial,
    FrobeniusCurveType, FrobeniusCurveTypeReport, FrobeniusDiscriminant,
    FrobeniusExtensionCountReport, FrobeniusExtensionCountSequenceReport,
    FrobeniusExtensionEnumerationComparisonReport, FrobeniusLocalZetaFunction,
    FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport, FrobeniusOrbit,
    FrobeniusTorsionMatrixError, FrobeniusTorsionMatrixReport, FrobeniusTrace,
    FrobeniusTraceCurveModel, GroupOrderReport, GroupOrderStrategy, HasseBoundReport,
    HasseGroupOrderStrategy, HasseInterval, HasseMultipleSearchReport, HasseMultipleSearchStep,
    IsogenyFrobeniusRelation, IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData,
    ModNMatrix2, NTorsionBasis, QuadraticTwistFrobeniusRelation, RelativeFrobenius,
    absolute_frobenius_on_exact_torsion, absolute_frobenius_orbit,
    absolute_frobenius_orbits_on_points, absolute_frobenius_power_point,
    compare_extension_count_with_enumeration, frobenius_matrix_on_n_torsion_basis,
    frobenius_twist_power, relative_frobenius_on_exact_torsion, relative_frobenius_orbit,
    relative_frobenius_orbits_on_points, relative_frobenius_point,
    verify_frobenius_characteristic_equation_at_point,
    verify_frobenius_characteristic_equation_exhaustive, verify_hasse_bound,
    verify_isogeny_frobenius_relation, verify_isogeny_graph_frobenius_relation,
};
pub use function_fields::{
    ShortWeierstrassFunction, ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
};
pub use invariants::HasJInvariant;
pub use isomorphisms::{
    CurveIsomorphism, CurveIsomorphismError, ShortWeierstrassIsomorphism,
    ShortWeierstrassQuadraticTwist, ShortWeierstrassTwist, TwistKind,
};
pub use short_weierstrass::{
    ExhaustivePointOrderReport, ExponentAccumulationReport, ExponentAccumulationStep,
    ExponentLowerBoundGroupOrderVerification, GroupExponentReport, GroupExponentStrategy,
    HasseIntervalPointOrderReport, PointOrderFromMultipleReport, PointOrderReductionStep,
    PointOrderReport, PointOrderStrategy, PointOrderStrategyKind, ShortWeierstrassCurve,
};
pub use torsion::{point_has_exact_order, points_of_exact_order};
pub use traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteAbelianGroupStructure,
    FiniteGroupCurveModel, GroupCurveModel, LiftXCoordinate, PointIndexSampler,
    RelativeFrobeniusCurveModel,
};
