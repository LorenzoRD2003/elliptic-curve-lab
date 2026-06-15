use crate::fields::traits::Field;
use crate::polynomials::PolynomialError;

use super::MultivariatePolynomial;

impl<F: Field> MultivariatePolynomial<F> {
    /// Evaluates a multivariate polynomial over a field at a concrete point.
    ///
    /// The point must provide exactly one coordinate for each variable in the
    /// polynomial's ambient arity.
    ///
    /// Each term is evaluated as
    ///
    /// `coefficient * x_0^e0 * x_1^e1 * ... * x_{n-1}^{e_{n-1}}`
    ///
    /// and the term contributions are then added together.
    pub fn evaluate(&self, point: &[F::Elem]) -> Result<F::Elem, PolynomialError> {
        if point.len() != self.arity() {
            return Err(PolynomialError::EvaluationPointArityMismatch {
                expected: self.arity(),
                actual: point.len(),
            });
        }

        let mut value = F::zero();

        for term in self.terms() {
            let mut monomial_value = F::one();

            for (coordinate, exponent) in point.iter().zip(&term.monomial.exponents) {
                let power = F::pow(coordinate, *exponent as u64);
                monomial_value = F::mul(&monomial_value, &power);
            }

            let contribution = F::mul(&term.coefficient, &monomial_value);
            value = F::add(&value, &contribution);
        }

        Ok(value)
    }
}
