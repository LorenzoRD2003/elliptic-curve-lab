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
    AnalyticTorsionPointApprox, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
    ComplexTorusPoint, EisensteinSumApprox, EllipticFunctionApproximation,
    EllipticFunctionTruncation, EvenDivisionPolynomialVanishingBranch,
    FundamentalParallelogramCoordinate, HasPoleDistance, LatticeIndexPoint, LatticeSumTruncation,
    TorusToCurveMapResult, TorusToCurveValues, TorusTorsionIndex, TorusTorsionPoint,
    TruncationConvergenceReport, UpperHalfPlanePoint, WeierstrassDifferentialEquationReport,
    WeierstrassDifferentialEquationStatus, WeierstrassPApprox, WeierstrassPDerivativeApprox,
    analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
    analytic_invariants_from_tau, analytic_j_invariant,
    compare_analytic_torsion_with_division_polynomial, compare_eisenstein_truncations,
    compare_primitive_analytic_torsion_with_division_polynomial, eisenstein as analytic_eisenstein,
    eisenstein_sum, elliptic_functions as analytic_elliptic_functions, errors as analytic_errors,
    evaluate_truncated_elliptic_function, explain as analytic_explain,
    fundamental_domain as analytic_fundamental_domain, g4_sum, g6_sum,
    invariants as analytic_invariants_module, lattice as analytic_lattice,
    map_fundamental_point_to_curve, map_primitive_torus_torsion_to_curve, map_torus_point_to_curve,
    map_torus_torsion_to_curve, modular_action as analytic_modular_action,
    periods as analytic_periods, primitive_torus_n_torsion_points,
    q_expansion as analytic_q_expansion, tolerance as analytic_tolerance,
    torsion as analytic_torsion, torus_n_torsion_points, torus_point as analytic_torus_point,
    truncation as analytic_truncation, upper_half_plane as analytic_upper_half_plane,
    verify_weierstrass_differential_equation, weierstrass_model as analytic_weierstrass_model,
    weierstrass_p, weierstrass_p_derivative,
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
