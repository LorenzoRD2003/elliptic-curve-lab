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

pub use elliptic_curves::{
    AffineCurveModel, AffinePoint, AnalyticCurveMembershipReport, AnalyticCurvePoint,
    AnalyticDivisionPolynomialComparisonCase, AnalyticDivisionPolynomialComparisonStatus,
    AnalyticEvenDivisionPolynomialReport, AnalyticInvariants, AnalyticOddDivisionPolynomialReport,
    AnalyticShortWeierstrassModel, AnalyticTorsionPointApprox, AnalyticWeierstrassCurve,
    ComplexAnalyticCurveLabReport, ComplexApproxComparison, ComplexDifferenceReport,
    ComplexLattice, ComplexTorusPoint, CurveError, CurveIsomorphism, CurveModel,
    EisensteinSeriesQExpansion, EisensteinSeriesQExpansionApprox, EisensteinSeriesWeight,
    EisensteinSumApprox, EllipticFunctionApproximation, EllipticFunctionTruncation,
    EnumerableCurveModel, EvenDivisionPolynomialVanishingBranch, FiniteAbelianGroupStructure,
    FiniteGroupCurveModel, FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
    FundamentalDomainReductionStep, FundamentalDomainReductionStepReason,
    FundamentalParallelogramCoordinate, GroupCurveModel, HasAnalyticLatticeContext,
    HasComplexApproxComparison, HasPoleDistance, JInvariantComparisonReport, JInvariantQExpansion,
    JInvariantQExpansionApprox, LatticeIndexPoint, LatticeSumTruncation, LiftXCoordinate,
    ModularInvarianceReport, ModularMatrix, ModularQExpansionApproximation,
    ModularQExpansionCoefficients, ModularQExpansionFamily, ModularQParameter, PointIndexSampler,
    QExpansionTruncation, ShortWeierstrassCurve, SpecialJKind, SpecialTauKind,
    TorusToCurveMapResult, TorusToCurveValues, TorusTorsionIndex, TorusTorsionPoint,
    TruncationConvergenceReport, UniformizationExperimentReport, UpperHalfPlanePoint,
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
    WeierstrassPApprox, WeierstrassPDerivativeApprox, WeierstrassZetaApprox, analytic_comparison,
    analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
    analytic_invariants_from_tau, analytic_j_invariant,
    compare_analytic_torsion_with_division_polynomial, compare_eisenstein_truncations,
    compare_j_from_eisenstein_and_q_expansion,
    compare_primitive_analytic_torsion_with_division_polynomial, eisenstein_sum,
    evaluate_truncated_elliptic_function, g4_sum, g6_sum, is_in_standard_fundamental_domain,
    map_fundamental_point_to_curve, map_primitive_torus_torsion_to_curve, map_torus_point_to_curve,
    map_torus_torsion_to_curve, point_has_exact_order, points_of_exact_order,
    primitive_torus_n_torsion_points, reduce_tau_to_standard_fundamental_domain,
    torus_n_torsion_points, verify_j_modular_invariance, verify_weierstrass_differential_equation,
    weierstrass_p, weierstrass_p_derivative, weierstrass_zeta,
};
pub use fields::{
    ApproxComparisonReport, ComplexApprox, EnumerableFiniteField, ExtensionField,
    ExtensionFieldElement, ExtensionFieldSpec, Field, FieldError, FiniteField,
    FiniteFieldDescriptor, Fp, FpElem, PolynomialFieldElement, PolynomialModulus, Q, SqrtField,
    addition_table, describe_complex, describe_complex_polynomial_modulus_as_field_modulus,
    describe_extension_field, describe_extension_field_element,
    describe_prime_polynomial_field_element, describe_prime_polynomial_modulus,
    describe_prime_polynomial_modulus_as_field_modulus, describe_rational, explain_add,
    explain_complex_polynomial_modulus_irreducibility, explain_complex_square_root,
    explain_extension_field_add, explain_extension_field_inverse, explain_extension_field_mul,
    explain_extension_field_reduction, explain_inverse, explain_mul,
    explain_prime_field_square_root, explain_prime_polynomial_field_add,
    explain_prime_polynomial_field_inverse, explain_prime_polynomial_field_mul,
    explain_prime_polynomial_field_reduction, explain_prime_polynomial_modulus_irreducibility,
    explain_prime_polynomial_storage, explain_rational_add, explain_rational_div,
    explain_rational_inverse, explain_rational_mul, explain_rational_square_root, format_complex,
    format_complex_polynomial, format_extension_field, format_extension_field_element,
    format_fp_elem, format_prime_field, format_prime_polynomial,
    format_prime_polynomial_field_element, format_prime_polynomial_modulus, format_rational,
    format_rational_field, inverses_table, multiplication_table,
};
pub use isogenies::{
    ComposedIsogeny, DualVeluIsogeny, Isogeny, IsogenyError, IsogenyKernel, IsomorphismIsogeny,
    ScalarMultiplicationIsogeny, VeluIsogeny, VerifiableIsogeny, maps_equal_exhaustively,
    verify_left_dual_relation, verify_right_dual_relation,
};
pub use numerics::{
    ApproxTolerance, bernoulli_number, sigma_power_sum_factorized, sigma_power_sum_naive,
    sigma_power_sums_up_to,
};
pub use polynomials::{
    DensePolynomial, IrreducibilityBackend, IrreducibilityStatus, PolynomialError,
    ReducibilityReason, SparsePolynomial, VisualizablePolynomial, describe_dense_polynomial,
    describe_multivariate_polynomial, describe_sparse_polynomial, explain_dense_storage,
    explain_multivariate_storage, explain_sparse_storage, format_dense_polynomial, format_monomial,
    format_multivariate_polynomial, format_sparse_polynomial, irreducibility_status,
    is_irreducible,
};
pub use visualization::{
    DivisionPolynomialKind, DivisionPolynomialSummary, IsogenyGraphSummary,
    VolcanoHeuristicSummary, describe_analytic_curve_membership,
    describe_analytic_division_polynomial_comparison,
    describe_analytic_even_division_polynomial_report, describe_analytic_invariants,
    describe_analytic_odd_division_polynomial_report, describe_analytic_torsion_point_approx,
    describe_complex_lattice, describe_composition, describe_curve, describe_dual_isogeny,
    describe_eisenstein_sum, describe_fundamental_domain_reduction_report,
    describe_fundamental_domain_reduction_step, describe_group_structure, describe_isogeny,
    describe_isomorphism, describe_j_invariant_comparison, describe_membership,
    describe_modular_invariance_report, describe_modular_matrix, describe_order_distribution,
    describe_point, describe_point_order, describe_q_parameter, describe_scalar_mul,
    describe_scalar_multiplication_isogeny, describe_torus_to_curve_map,
    describe_truncation_convergence, describe_weierstrass_differential_equation,
    describe_weierstrass_p_approx, describe_weierstrass_p_derivative_approx,
    division_polynomial_summary, explain_add as explain_curve_add, explain_division_polynomial,
    explain_dual_relation, explain_isogeny_graph, explain_point_order, explain_quadratic_twist,
    explain_short_weierstrass_scaling, explain_torsion_via_division_polynomial,
    explain_velu_codomain, explain_velu_evaluation, explain_volcano_like_layers,
    format_adjacency_list, format_analytic_cubic_model, format_curve, format_isogeny,
    format_isomorphism, format_point, format_point_compact, format_short_weierstrass_over_complex,
    list_points, summarize_curve_comparison, summarize_dual_verification,
    summarize_group_structure, summarize_kernel, summarize_order_distribution,
};
pub use visualization::{Visualizable, VisualizableField};
