use crate::fields::traits::Field;
use crate::polynomials::PolynomialError;

use super::SparsePolynomial;

impl<F: Field> SparsePolynomial<F> {
    /// Evaluates a sparse univariate polynomial over a field at a point of the
    /// same field.
    ///
    /// Because the representation stores only non-zero terms, the implementation
    /// evaluates each stored term as
    ///
    /// `coefficient * point^degree`
    ///
    /// and sums the results.
    pub fn evaluate(&self, point: &F::Elem) -> Result<F::Elem, PolynomialError> {
        let mut value = F::zero();

        for term in self.terms() {
            let power = F::pow(point, term.degree as u64);
            let contribution = F::mul(&term.coefficient, &power);
            value = F::add(&value, &contribution);
        }

        Ok(value)
    }
}
