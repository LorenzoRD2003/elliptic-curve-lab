//! Elliptic-curve scaffolding.

pub mod affine;
pub mod division_polynomials;
pub mod error;
pub mod invariants;
pub mod isomorphisms;
pub mod short_weierstrass;
pub mod torsion;
pub mod traits;

pub use affine::AffinePoint;
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
