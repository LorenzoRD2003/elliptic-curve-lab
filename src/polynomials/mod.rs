//! Polynomial data structures and algorithms.

pub mod dense;
pub mod error;
pub mod irreducibility;
pub mod multivariate;
pub mod sparse;
pub mod traits;

pub use dense::DensePolynomial;
pub use error::PolynomialError;
pub use multivariate::MultivariatePolynomial;
pub use sparse::SparsePolynomial;
