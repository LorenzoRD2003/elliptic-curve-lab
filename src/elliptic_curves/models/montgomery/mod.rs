//! Montgomery curves.
//!
//! The target is the affine model
//!
//! `B y^2 = x^3 + A x^2 + x`
//!
//! over fields of characteristic different from `2`.

mod api;
mod display;
mod group_law;
mod invariants;
mod membership;
mod model_traits;
mod reduction;
mod type_definition;

#[cfg(test)]
mod tests;

pub use type_definition::MontgomeryCurve;
