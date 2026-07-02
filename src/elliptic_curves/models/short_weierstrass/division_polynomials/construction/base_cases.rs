use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::traits::*;
use crate::polynomials::DensePolynomial;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Returns the base short-Weierstrass division polynomials `ψ_0` through `ψ_4`.
    ///
    /// This helper currently covers only the explicit low-degree formulas:
    ///
    /// - `ψ_0 = 0`
    /// - `ψ_1 = 1`
    /// - `ψ_2 = 2y`
    /// - `ψ_3 = 3x^4 + 6ax^2 + 12bx - a^2`
    /// - `ψ_4 = 4y (x^6 + 5ax^4 + 20bx^3 - 5a^2x^2 - 4abx - 8b^2 - a^3)`
    ///
    /// For larger indices, the current division-polynomial scaffold reports
    /// [`DivisionPolynomialError::UnsupportedIndex`].
    ///
    /// Complexity: `Θ(1)` field operations and bounded-size coefficient vectors.
    pub(crate) fn base_division_polynomial(
        &self,
        n: usize,
    ) -> Result<DivisionPolynomialForm<F>, DivisionPolynomialError> {
        fn scaled<F: Field>(coefficient: &F::Elem, scalar: i64) -> F::Elem {
            F::mul(coefficient, &F::from_i64(scalar))
        }

        let a = self.a();
        let b = self.b();
        let a2 = F::square(a);
        let a3 = F::mul(&a2, a);
        let b2 = F::square(b);
        let ab = F::mul(a, b);

        match n {
            0 => Ok(DivisionPolynomialForm::x_polynomial(DensePolynomial::new(
                Vec::new(),
            ))),
            1 => Ok(DivisionPolynomialForm::x_polynomial(
                DensePolynomial::constant(F::one()),
            )),
            2 => Ok(DivisionPolynomialForm::y_times_x_polynomial(
                DensePolynomial::constant(F::from_i64(2)),
            )),
            3 => Ok(DivisionPolynomialForm::x_polynomial(DensePolynomial::new(
                vec![
                    F::neg(&a2),
                    scaled::<F>(b, 12),
                    scaled::<F>(a, 6),
                    F::zero(),
                    F::from_i64(3),
                ],
            ))),
            4 => {
                let constant_term = F::add(&scaled::<F>(&b2, -32), &scaled::<F>(&a3, -4));

                Ok(DivisionPolynomialForm::y_times_x_polynomial(
                    DensePolynomial::new(vec![
                        constant_term,
                        scaled::<F>(&ab, -16),
                        scaled::<F>(&a2, -20),
                        scaled::<F>(b, 80),
                        scaled::<F>(a, 20),
                        F::zero(),
                        F::from_i64(4),
                    ]),
                ))
            }
            _ => Err(DivisionPolynomialError::UnsupportedIndex { n }),
        }
    }
}
