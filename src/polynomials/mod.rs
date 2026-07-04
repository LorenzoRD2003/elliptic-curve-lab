//! Polynomial data structures and algorithms.

pub mod dense;
pub mod error;
pub(crate) mod integer;
pub mod irreducibility;
pub mod multivariate;
pub(crate) mod rational_normalization;
pub mod sparse;
pub mod traits;

pub use dense::DensePolynomial;
pub use error::PolynomialError;
pub(crate) use integer::IntegerPolynomial;
pub use multivariate::MultivariatePolynomial;
pub(crate) use rational_normalization::primitive_integer_polynomial;
pub use sparse::SparsePolynomial;
