pub mod dense;
pub mod division;
pub mod evaluation;
pub mod gcd;
pub mod interpolation;
pub mod irreducibility;
pub mod multivariate;
pub mod sparse;
pub mod traits;

pub use dense::{describe_dense_polynomial, explain_dense_storage, format_dense_polynomial};
pub use division::explain_dense_division;
pub use evaluation::{
    explain_evaluate_dense, explain_evaluate_multivariate, explain_evaluate_sparse,
};
pub use gcd::explain_dense_gcd;
pub use interpolation::explain_lagrange_interpolation;
pub use irreducibility::{describe_irreducibility_status, explain_dense_irreducibility};
pub use multivariate::{
    describe_multivariate_polynomial, explain_multivariate_storage, format_monomial,
    format_multivariate_polynomial,
};
pub use sparse::{describe_sparse_polynomial, explain_sparse_storage, format_sparse_polynomial};
pub use traits::VisualizablePolynomial;
