use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::division_polynomials::DivisionPolynomialError,
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassCurve<F> {
    pub(crate) fn division_polynomial_uses_y_factor(
        &self,
        n: usize,
    ) -> Result<bool, DivisionPolynomialError> {
        self.ensure_division_index_nonzero(n)?;
        Ok(n.is_multiple_of(2))
    }

    pub(crate) fn ensure_division_index_nonzero(
        &self,
        n: usize,
    ) -> Result<(), DivisionPolynomialError> {
        if n == 0 {
            Err(DivisionPolynomialError::ZeroIndex)
        } else {
            Ok(())
        }
    }

    pub(crate) fn ensure_odd_division_index(
        &self,
        n: usize,
    ) -> Result<(), DivisionPolynomialError> {
        self.ensure_division_index_nonzero(n)?;
        if n.is_multiple_of(2) {
            Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n })
        } else {
            Ok(())
        }
    }

    pub(crate) fn ensure_even_division_index(
        &self,
        n: usize,
    ) -> Result<(), DivisionPolynomialError> {
        self.ensure_division_index_nonzero(n)?;
        if n.is_multiple_of(2) {
            Ok(())
        } else {
            Err(DivisionPolynomialError::UnsupportedIndex { n })
        }
    }
}
