//! Polynomial data structures and algorithms.

pub mod dense;
pub mod error;
pub mod evaluation;
pub mod interpolation;
pub mod irreducibility;
pub mod multivariate;
pub mod sparse;
pub mod traits;

pub use dense::DensePolynomial;
pub use error::PolynomialError;
pub use irreducibility::{
    IrreducibilityBackend, IrreducibilityStatus, ReducibilityReason, irreducibility_status,
    is_irreducible,
};
pub use multivariate::{Monomial, MultivariatePolynomial, MultivariateTerm};
pub use sparse::{SparsePolynomial, SparsePolynomialTerm};
pub use traits::UnivariatePolynomial;
