use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::DivisionPolynomialError,
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

impl<F: EnumerableFiniteField + SqrtField> ShortWeierstrassCurve<F> {
    /// Returns the rational affine points found from the division-polynomial
    /// vanishing condition `ψ_n(P) = 0`.
    pub fn torsion_candidates_from_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
        if self.division_polynomial_uses_y_factor(n)? {
            self.torsion_candidates_from_even_division_polynomial(n)
        } else {
            self.torsion_candidates_from_odd_division_polynomial(n)
        }
    }

    pub(crate) fn torsion_candidates_from_odd_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
        self.ensure_odd_division_index(n)?;

        let mut points = Vec::new();
        for x in F::elements() {
            if !self.x_criterion_vanishes(n, &x)? {
                continue;
            }

            self.for_each_rational_affine_point_with_x(&x, |point| points.push(point.clone()));
        }

        Ok(points)
    }

    pub(crate) fn torsion_candidates_from_even_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<AffinePoint<F>>, DivisionPolynomialError> {
        self.ensure_even_division_index(n)?;

        let mut points = Vec::new();
        for x in F::elements() {
            let even_factor_zero = self.x_criterion_vanishes(n, &x)?;
            self.for_each_rational_affine_point_with_x(&x, |point| {
                if self.point_has_zero_y(point) || even_factor_zero {
                    points.push(point.clone());
                }
            });
        }

        Ok(points)
    }
}
