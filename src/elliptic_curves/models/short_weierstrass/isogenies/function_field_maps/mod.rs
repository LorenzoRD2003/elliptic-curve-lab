//! Pullback maps on short-Weierstrass function fields.
//!
//! For a morphism of short-Weierstrass curves `φ: E -> E'`,
//! the induced map on function fields goes in the opposite direction:
//!
//! `φ^* : F(E') -> F(E)`.
//!
//! In the current implementation, such a pullback is represented by the images
//! of the codomain coordinate functions `x'` and `y'` inside the domain function field `F(E)`.

mod map;
mod report;

#[cfg(test)]
mod tests;

pub use map::ShortWeierstrassFunctionFieldMap;
pub use report::{DifferentialPullbackReport, IsogenySeparabilityKind};
