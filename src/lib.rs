//! Foundational scaffolding for mathematical and cryptographic algorithms.
//!
//! The crate intentionally starts with small, documented interfaces and
//! lightweight placeholder implementations so the core abstractions can evolve
//! with tests before the heavy algebraic algorithms arrive.

// pub mod algorithms;
pub mod elliptic_curves;
pub mod fields;
pub mod isogenies;
pub mod numerics;
pub mod polynomials;
pub mod visualization;

#[cfg(any(test, feature = "test-support"))]
pub(crate) mod proptest_support;

pub use elliptic_curves::{
    AbelJacobiConfig, AbsoluteFrobenius, AffineCurveModel, AffinePoint, AnalyticCurveError,
    AnalyticCurvePoint, AnalyticInvariants, AnalyticShortWeierstrassModel,
    AnalyticWeierstrassCurve, CompleteEllipticIntegralKApprox, ComplexLattice,
    ComplexModuloLatticeComparison, ComplexTorusPoint, CubicRootConfiguration, CurveError,
    CurveIsomorphism, CurveModel, DiscriminantSign, EisensteinSeriesQExpansion,
    EisensteinSeriesWeight, EllipticFunctionTruncation, EndomorphismRingCandidateSet,
    EndomorphismRingLocalView, EndomorphismRingReport, EnumerableCurveModel,
    FiniteAbelianGroupStructure, FiniteGroupCurveModel, FrobeniusCharacteristicEquationCheck,
    FrobeniusCharacteristicEquationExhaustiveReport, FrobeniusCharacteristicPolynomial,
    FrobeniusCurveType, FrobeniusCurveTypeReport, FrobeniusDiscriminant,
    FrobeniusExtensionCountReport, FrobeniusExtensionCountSequenceReport,
    FrobeniusExtensionEnumerationComparisonReport, FrobeniusLocalZetaFunction,
    FrobeniusOnExactTorsionPoint, FrobeniusOnExactTorsionReport, FrobeniusOrbit,
    FrobeniusTorsionMatrixError, FrobeniusTorsionMatrixReport, FrobeniusTrace,
    FrobeniusTraceCurveModel, FundamentalParallelogramCoordinate, GroupCurveModel,
    HasseBoundReport, ImaginaryQuadraticOrder, ImaginaryQuadraticOrderError,
    IsogenyFrobeniusRelation, IsogenyGraphFrobeniusReport, IsogenyGraphNodeFrobeniusData,
    LatticeIndexPoint, LatticeSumTruncation, LegendreParameter, LegendreReduction, LiftXCoordinate,
    ModNMatrix2, ModularMatrix, ModularQParameter, NTorsionBasis, PeriodRecoveryConfig,
    PointIndexSampler, QExpansionTruncation, QuadraticDiscriminant,
    QuadraticDiscriminantFactorization, QuadraticDiscriminantFactorizationError,
    QuadraticDiscriminantMod4, QuadraticOrderCoverRelation, QuadraticOrderIndexError,
    QuadraticTwistFrobeniusRelation, RecoveredPeriodBasis, RelativeFrobenius,
    RelativeFrobeniusCurveModel, ShortWeierstrassCurve, TorusTorsionIndex, TorusTorsionPoint,
    UpperHalfPlanePoint, VolcanoEndomorphismLevelCandidate, WeierstrassCubicRoots,
    absolute_frobenius_on_exact_torsion, absolute_frobenius_orbit,
    absolute_frobenius_orbits_on_points, absolute_frobenius_power_point, analytic_discriminant,
    analytic_g2, analytic_g3, analytic_invariants, analytic_invariants_from_tau,
    analytic_j_invariant, approximate_abel_jacobi_integral,
    compare_analytic_torsion_with_division_polynomial, compare_extension_count_with_enumeration,
    compare_primitive_analytic_torsion_with_division_polynomial,
    complementary_complete_elliptic_integral_k_from_lambda,
    complementary_complete_elliptic_integral_k_from_m, complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_m, eisenstein_sum, frobenius_matrix_on_n_torsion_basis,
    frobenius_twist_power, g4_sum, g6_sum, map_fundamental_point_to_curve,
    map_primitive_torus_torsion_to_curve, map_torus_point_to_curve, map_torus_torsion_to_curve,
    point_has_exact_order, points_of_exact_order, primitive_torus_n_torsion_points,
    recover_canonical_tau_from_curve, recover_period_basis,
    recover_period_basis_from_legendre_reduction, recover_tau_from_curve,
    recover_torus_point_from_curve_point, recover_torus_point_from_curve_point_with_periods,
    recover_weierstrass_cubic_roots, recover_weierstrass_cubic_roots_from_invariants,
    reduce_tau_to_standard_fundamental_domain, relative_frobenius_on_exact_torsion,
    relative_frobenius_orbit, relative_frobenius_orbits_on_points, relative_frobenius_point,
    torus_n_torsion_points, verify_frobenius_characteristic_equation_at_point,
    verify_frobenius_characteristic_equation_exhaustive, verify_hasse_bound,
    verify_isogeny_frobenius_relation, verify_isogeny_graph_frobenius_relation, weierstrass_p,
    weierstrass_p_derivative, weierstrass_zeta,
};
pub use fields::{
    ComplexApprox, EnumerableFiniteField, ExtensionField, ExtensionFieldElement,
    ExtensionFieldSpec, Field, FieldError, FiniteField, Fp, FpElem, PolynomialFieldElement,
    PolynomialModulus, Q, SqrtField,
};
pub use isogenies::{
    ComposedIsogeny, DualVeluIsogeny, EndomorphismVolcanoReport, Isogeny,
    IsogenyEdgeEndomorphismRelation, IsogenyEdgeEndomorphismReport, IsogenyError,
    IsogenyGraphEndomorphismEdgeReport, IsogenyGraphEndomorphismNodeReport,
    IsogenyGraphEndomorphismReport, IsogenyKernel, IsomorphismIsogeny, ScalarMultiplicationIsogeny,
    VeluIsogeny, VerifiableIsogeny, VolcanoHeuristicComparison, maps_equal_exhaustively,
    verify_left_dual_relation, verify_right_dual_relation,
};
pub use numerics::{ApproxTolerance, PositivePrimeError, valuation_biguint};
pub use polynomials::{
    DensePolynomial, IrreducibilityBackend, IrreducibilityStatus, PolynomialError,
    ReducibilityReason, SparsePolynomial, irreducibility_status, is_irreducible,
};
