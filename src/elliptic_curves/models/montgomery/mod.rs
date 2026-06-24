//! Montgomery curves.
//!
//! The target is the affine model
//!
//! `B y^2 = x^3 + A x^2 + x`
//!
//! over fields of characteristic different from `2`.

mod api;
mod differential_arithmetic;
mod display;
mod group_law;
mod invariants;
mod membership;
mod model_traits;
mod normalization;
mod reduction;
mod type_definition;
mod x_coordinates;

#[cfg(test)]
mod tests;

pub use differential_arithmetic::MontgomeryDifferentialArithmeticError;
pub use normalization::{MontgomeryNormalizationError, NormalizedMontgomeryCurve};
pub use type_definition::MontgomeryCurve;
pub use x_coordinates::MontgomeryXzPoint;
