use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Returns the short-Weierstrass division polynomial in its honest public shape.
    pub fn division_polynomial(
        &self,
        n: usize,
    ) -> Result<DivisionPolynomialForm<F>, DivisionPolynomialError> {
        match n {
            0..=4 => self.base_division_polynomial(n),
            _ if n.is_multiple_of(2) => Ok(DivisionPolynomialForm::YTimes(
                self.even_division_polynomial_factor(n)?,
            )),
            _ => Ok(DivisionPolynomialForm::InX(
                self.odd_division_polynomial(n)?,
            )),
        }
    }
}
