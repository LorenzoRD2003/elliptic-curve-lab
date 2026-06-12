//! Pullback maps on short-Weierstrass function fields.
//!
//! For a morphism of short-Weierstrass curves `phi : E -> E'`,
//! the induced map on function fields goes in the opposite direction:
//!
//! `phi^* : F(E') -> F(E)`.
//!
//! In the current implementation, such a pullback is represented by the images
//! of the codomain coordinate functions `x'` and `y'` inside the domain function field `F(E)`.

mod map;
mod report;

#[cfg(test)]
mod tests;

pub use map::ShortWeierstrassFunctionFieldMap;
pub(crate) use report::DifferentialPullbackReportParts;
pub use report::{DifferentialPullbackReport, IsogenySeparabilityKind};
