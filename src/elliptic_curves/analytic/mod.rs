//! Educational scaffolding for the future complex-analytic elliptic-curve
//! milestone.

pub mod eisenstein;
pub mod elliptic_functions;
pub mod errors;
pub mod explain;
pub mod fundamental_domain;
pub mod invariants;
pub mod lattice;
pub mod modular_action;
pub mod periods;
pub mod q_expansion;
pub mod reports;
pub mod torsion;
pub mod torus_point;
pub mod truncation;
pub mod upper_half_plane;
pub mod weierstrass_model;

pub use crate::numerics::ApproxTolerance;
pub use crate::numerics::tolerance;
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
pub use invariants::{
    AnalyticInvariants, analytic_discriminant, analytic_g2, analytic_g3, analytic_invariants,
    analytic_invariants_from_tau, analytic_j_invariant,
};
pub use lattice::{
    ComplexLattice, ComplexTorusPoint, FundamentalParallelogramCoordinate, LatticeIndexPoint,
};
pub use reports::{
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
    verify_weierstrass_differential_equation,
};
pub use torus_point::{
    TorusToCurveMapResult, TorusToCurveValues, map_fundamental_point_to_curve,
    map_torus_point_to_curve,
};
pub use truncation::LatticeSumTruncation;
pub use upper_half_plane::UpperHalfPlanePoint;
pub use weierstrass_model::{
    AnalyticCurveMembershipReport, AnalyticCurvePoint, AnalyticWeierstrassCurve,
};
