use crate::fields::traits::Field;
use crate::polynomials::PolynomialError;

use super::DensePolynomial;

impl<F: Field> DensePolynomial<F> {
    /// Evaluates a dense univariate polynomial over a field at a point of the same
    /// field.
    ///
    /// The implementation uses Horner's rule on the coefficient vector stored in
    /// ascending degree order. This gives a clear and efficient baseline while
    /// staying easy to explain.
    ///
    /// If the polynomial is represented by an empty coefficient vector, the result
    /// is the additive identity of the field.
    pub fn evaluate(&self, point: &F::Elem) -> Result<F::Elem, PolynomialError> {
        let mut accumulator = F::zero();

        for coefficient in self.coefficients().iter().rev() {
            accumulator = F::add(&F::mul(&accumulator, point), coefficient);
        }

        Ok(accumulator)
    }
}
