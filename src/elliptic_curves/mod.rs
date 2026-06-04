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
    AnalyticCurveError, AnalyticCurveMembershipReport, AnalyticCurvePoint, AnalyticInvariants,
    AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice, ComplexTorusPoint,
    EisensteinSumApprox, FundamentalParallelogramCoordinate, LatticeIndexPoint,
    LatticeSumTruncation, TruncationConvergenceReport, UpperHalfPlanePoint, analytic_discriminant,
    analytic_g2, analytic_g3, analytic_invariants, analytic_invariants_from_tau,
    analytic_j_invariant, compare_eisenstein_truncations, eisenstein as analytic_eisenstein,
    eisenstein_sum, elliptic_functions as analytic_elliptic_functions, errors as analytic_errors,
    explain as analytic_explain, fundamental_domain as analytic_fundamental_domain, g4_sum, g6_sum,
    invariants as analytic_invariants_module, lattice as analytic_lattice,
    modular_action as analytic_modular_action, periods as analytic_periods,
    q_expansion as analytic_q_expansion, reports as analytic_reports,
    tolerance as analytic_tolerance, torsion as analytic_torsion,
    torus_point as analytic_torus_point, truncation as analytic_truncation,
    upper_half_plane as analytic_upper_half_plane, weierstrass_model as analytic_weierstrass_model,
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
