use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::DivisionPolynomialError,
    traits::{CurveModel, LiftXCoordinate},
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

impl<F: EnumerableFiniteField + SqrtField> ShortWeierstrassCurve<F> {
    /// Returns the rational `x`-coordinates in the base field that can correspond
    /// to rational affine points annihilated by the division polynomial `ψ_n`.
    pub fn rational_x_candidates_for_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<F::Elem>, DivisionPolynomialError> {
        if self.division_polynomial_uses_y_factor(n)? {
            self.rational_x_candidates_from_even_division_polynomial(n)
        } else {
            self.rational_roots_of_odd_division_polynomial(n)
        }
    }

    pub(crate) fn rational_roots_of_odd_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<F::Elem>, DivisionPolynomialError> {
        self.ensure_odd_division_index(n)?;

        let mut roots = Vec::new();
        for x in F::elements() {
            if self.x_criterion_vanishes(n, &x)? && self.point_from_x(x.clone())?.is_some() {
                roots.push(x);
            }
        }
        Ok(roots)
    }

    pub(crate) fn rational_x_candidates_from_even_division_polynomial(
        &self,
        n: usize,
    ) -> Result<Vec<F::Elem>, DivisionPolynomialError> {
        self.ensure_even_division_index(n)?;

        let mut xs = Vec::new();
        for x in F::elements() {
            let mut should_keep = false;
            self.for_each_rational_affine_point_with_x(&x, |point| {
                should_keep = self.point_satisfies_even_division_x_candidate_condition(n, point);
            });
            if should_keep {
                xs.push(x);
            }
        }

        Ok(xs)
    }

    pub(crate) fn x_criterion_vanishes(
        &self,
        n: usize,
        x: &F::Elem,
    ) -> Result<bool, DivisionPolynomialError> {
        self.evaluate_division_polynomial_x_criterion(n, x)
            .map(|value| F::is_zero(&value))
    }

    fn point_satisfies_even_division_x_candidate_condition(
        &self,
        n: usize,
        point: &AffinePoint<F>,
    ) -> bool {
        self.point_has_zero_y(point)
            || self
                .finite_point_x(point)
                .is_some_and(|x| self.x_criterion_vanishes(n, x).is_ok_and(|value| value))
    }

    pub(crate) fn point_has_zero_y(&self, point: &AffinePoint<F>) -> bool {
        match point {
            AffinePoint::Finite { y, .. } => F::is_zero(y),
            AffinePoint::Infinity => false,
        }
    }

    pub(crate) fn finite_point_x<'a>(&self, point: &'a AffinePoint<F>) -> Option<&'a F::Elem> {
        match point {
            AffinePoint::Finite { x, .. } => Some(x),
            AffinePoint::Infinity => None,
        }
    }

    pub(super) fn for_each_rational_affine_point_with_x(
        &self,
        x: &F::Elem,
        mut visit: impl FnMut(&AffinePoint<F>),
    ) {
        for y in F::elements() {
            let point = AffinePoint::<F>::new(x.clone(), y);
            if self.contains(&point) {
                visit(&point);
            }
        }
    }
}
