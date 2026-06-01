//! Educational scaffolding for explicit curve isomorphisms.
//!
//! For the short-Weierstrass milestone, the intended convention is the
//! coordinate change
//! `\phi_u : E -> E'`,
//! given on affine points by
//! `(x, y) -> (u^2 x, u^3 y)`.
//!
//! When
//! `E : y^2 = x^3 + ax + b`,
//! this convention sends the equation to another short-Weierstrass model
//! `E' : y^2 = x^3 + a'x + b'`
//! with
//! `a' = u^4 a`
//! and
//! `b' = u^6 b`.
//!
//! Future implementations in this module should treat that normalization as the
//! canonical meaning of a short-Weierstrass scaling isomorphism.
mod error;
mod short_weierstrass;

pub use error::CurveIsomorphismError;
pub use short_weierstrass::{
    ShortWeierstrassIsomorphism, ShortWeierstrassQuadraticTwist, ShortWeierstrassTwist, TwistKind,
};
