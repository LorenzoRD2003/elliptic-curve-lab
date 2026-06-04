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
pub mod upper_half_plane;
pub mod weierstrass_model;

pub use crate::numerics::ApproxTolerance;
pub use crate::numerics::tolerance;
pub use errors::AnalyticCurveError;
pub use lattice::{
    ComplexLattice, ComplexTorusPoint, FundamentalParallelogramCoordinate, LatticeIndexPoint,
};
pub use upper_half_plane::UpperHalfPlanePoint;
