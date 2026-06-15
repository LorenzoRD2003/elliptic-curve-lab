//! Multivariate polynomials over a field.

mod arithmetic;
mod evaluation;
mod monomial;
mod tests;
mod type_definition;

pub use monomial::Monomial;
pub use type_definition::{MultivariatePolynomial, MultivariateTerm};
