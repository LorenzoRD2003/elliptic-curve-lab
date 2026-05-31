//! Elliptic-curve scaffolding.

pub mod affine;
pub mod error;
pub mod short_weierstrass;
pub mod traits;

pub use affine::AffinePoint;
pub use error::CurveError;
pub use short_weierstrass::ShortWeierstrassCurve;
pub use traits::{
    AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
    LiftXCoordinate, PointIndexSampler,
};
