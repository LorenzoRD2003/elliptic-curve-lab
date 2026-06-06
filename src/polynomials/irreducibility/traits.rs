use crate::fields::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};

use crate::polynomials::irreducibility::IrreducibilityStatus;

/// Backend capability marker for dense irreducibility classification.
///
/// The public irreducibility API is generic over this trait instead of naming
/// a specific algorithm in the function names. That lets the crate evolve from
/// the current educational baselines to stronger backend-specific algorithms
/// without forcing public renames.
pub trait IrreducibilityBackend: Field + Sized {
    /// Returns a structured irreducibility classification for dense
    /// polynomials over this base-field backend.
    fn irreducibility_status_impl(
        polynomial: &DensePolynomial<Self>,
    ) -> Result<IrreducibilityStatus<Self>, PolynomialError>;
}
