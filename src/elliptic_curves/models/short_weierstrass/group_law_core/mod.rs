//! Shared affine short-Weierstrass formulas across coordinate backends.
//!
//! The current short-Weierstrass model uses the same secant/tangent formulas
//! over two ambient coordinate stories:
//!
//! - ordinary affine points over the base field
//! - generic/function-field points over `F(E)`
//!
//! This module keeps the shared formula logic behind a small internal point
//! shape, a backend-operations trait, and one runner that owns the curve
//! coefficient `a` together with the chosen coordinate operations.

mod formulas;
mod ops;
mod point;
mod runner;

#[cfg(test)]
mod tests;

pub(crate) use ops::ShortWeierstrassFormulaOps;
pub(crate) use point::ShortWeierstrassFormulaPoint;
pub(crate) use runner::ShortWeierstrassFormulaRunner;
