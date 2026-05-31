//! Polynomial data structures and algorithms.

pub mod dense;
pub mod error;
pub mod evaluation;
pub mod interpolation;
pub mod multivariate;
pub mod sparse;
pub mod traits;

pub use crate::visualization::polynomials::{
    VisualizablePolynomial, describe_dense_polynomial, describe_multivariate_polynomial,
    describe_sparse_polynomial, explain_dense_division, explain_dense_gcd, explain_dense_storage,
    explain_evaluate_dense, explain_evaluate_multivariate, explain_evaluate_sparse,
    explain_lagrange_interpolation, explain_multivariate_storage, explain_sparse_storage,
    format_dense_polynomial, format_monomial, format_multivariate_polynomial,
    format_sparse_polynomial,
};
pub use dense::DensePolynomial;
pub use error::PolynomialError;
pub use multivariate::{Monomial, MultivariatePolynomial, MultivariateTerm};
pub use sparse::{SparsePolynomial, SparsePolynomialTerm};
pub use traits::UnivariatePolynomial;
