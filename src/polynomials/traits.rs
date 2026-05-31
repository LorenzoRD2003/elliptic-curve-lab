use crate::fields::Field;
use crate::polynomials::PolynomialError;

/// Shared algebraic interface for univariate polynomials over a field.
///
/// The goal of this trait is to capture the common univariate-polynomial
/// operations that make sense across multiple concrete representations such as
/// dense and sparse storage.
///
/// This trait is intentionally smaller than any one concrete implementation:
///
/// - it does not expose representation-specific accessors like
///   `coefficients()` or `terms()`
/// - it does not require Euclidean division, because that is currently only
///   implemented for dense polynomials
/// - it focuses on the common algebraic surface
pub trait UnivariatePolynomial<F: Field>: Sized {
    /// Builds the constant polynomial with the given value.
    fn constant(value: F::Elem) -> Self;

    /// Returns the degree of the polynomial when it is non-zero.
    fn degree(&self) -> Option<usize>;

    /// Returns the leading coefficient when the polynomial is non-zero.
    fn leading_coefficient(&self) -> Option<&F::Elem>;

    /// Returns the constant term if it is explicitly present.
    fn constant_term(&self) -> Option<&F::Elem>;

    /// Returns whether the polynomial is the zero polynomial.
    fn is_zero(&self) -> bool;

    /// Returns whether the polynomial is monic.
    ///
    /// By default, a polynomial is considered monic when it has a leading
    /// coefficient and that coefficient equals `1` in the base field.
    fn is_monic(&self) -> bool {
        self.leading_coefficient()
            .is_some_and(|leading| F::eq(leading, &F::one()))
    }

    /// Adds two univariate polynomials.
    fn add(&self, rhs: &Self) -> Self;

    /// Negates every coefficient of the polynomial.
    fn neg(&self) -> Self;

    /// Subtracts two univariate polynomials.
    fn sub(&self, rhs: &Self) -> Self;

    /// Multiplies every coefficient by the same field element.
    fn scale(&self, scalar: &F::Elem) -> Self;

    /// Multiplies two univariate polynomials.
    fn mul(&self, rhs: &Self) -> Self;

    /// Returns the monic normalization of the polynomial.
    ///
    /// The default implementation scales the polynomial by the inverse of its
    /// leading coefficient. The zero polynomial has no monic normalization.
    fn make_monic(&self) -> Result<Self, PolynomialError> {
        let Some(leading) = self.leading_coefficient() else {
            return Err(PolynomialError::ZeroPolynomialHasNoMonicNormalization);
        };

        let inverse =
            F::inverse(leading).map_err(|_| PolynomialError::NonInvertibleLeadingCoefficient)?;
        Ok(self.scale(&inverse))
    }
}
