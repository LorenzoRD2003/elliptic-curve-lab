//! Sparse univariate polynomials over a field.

mod arithmetic;
mod conversions;
mod evaluation;
mod tests;
mod trait_impls;
mod type_definition;

pub use type_definition::{SparsePolynomial, SparsePolynomialTerm};
