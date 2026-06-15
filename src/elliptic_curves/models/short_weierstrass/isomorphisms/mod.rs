//! Educational scaffolding for explicit short-Weierstrass curve isomorphisms.
//!
//! For the short-Weierstrass isomorphism layer, the intended convention is the
//! coordinate change  `ϕ_u : E -> E'`, given by `(x, y) -> (u^2 x, u^3 y)`.
//!
//! When `E : y^2 = x^3 + ax + b`, this convention sends the equation to another
//! short-Weierstrass model `E' : y^2 = x^3 + a'x + b'` with `a' = u^4 a` and
//! `b' = u^6 b`.
//!
//! Future implementations in this module should treat that normalization as the
//! canonical meaning of a short-Weierstrass scaling isomorphism.
mod error;
mod scaling;
mod twist;

#[cfg(test)]
mod tests;

pub use error::CurveIsomorphismError;
pub use scaling::ShortWeierstrassIsomorphism;
pub use twist::{ShortWeierstrassQuadraticTwist, ShortWeierstrassTwist, TwistKind};
