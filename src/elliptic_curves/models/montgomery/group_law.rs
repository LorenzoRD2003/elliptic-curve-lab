use crate::elliptic_curves::{
    AffinePoint, CurveError, MontgomeryCurve,
    traits::{CurveModel, GroupCurveModel},
};
use crate::fields::traits::Field;

impl<F: Field> MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    // TODO: the long-term executable core for Montgomery arithmetic should be
    // a dedicated projective-coordinate group law rather than this staged
    // affine implementation.

    fn affine_negated_y(&self, y: &F::Elem) -> F::Elem {
        F::neg(y)
    }

    /// Returns the affine tangent slope numerator
    ///
    /// `3x^2 + 2Ax + 1`
    ///
    /// for the Montgomery model `B y^2 = x^3 + A x^2 + x`.
    fn doubling_slope_numerator(&self, x: &F::Elem) -> F::Elem {
        F::add(
            &F::add(
                &F::mul(&F::from_i64(3), &F::square(x)),
                &F::mul(&F::from_i64(2), &F::mul(self.a(), x)),
            ),
            &F::one(),
        )
    }

    /// Adds two already-validated finite affine points with distinct
    /// `x`-coordinates.
    ///
    /// For the Montgomery model `B y^2 = x^3 + A x^2 + x`, the secant slope is
    ///
    /// `λ = (y2 - y1) / (x2 - x1)`
    ///
    /// and the resulting sum has coordinates
    ///
    /// `x3 = B λ^2 - A - x1 - x2`, `y3 = λ (x1 - x3) - y1`.
    fn add_distinct_finite_points(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> Result<AffinePoint<F>, CurveError> {
        let slope = F::div(&F::sub(right_y, left_y), &F::sub(right_x, left_x))
            .expect("distinct x-coordinates give a nonzero secant denominator");
        let x3 = F::sub(
            &F::sub(
                &F::sub(&F::mul(self.b(), &F::square(&slope)), self.a()),
                left_x,
            ),
            right_x,
        );
        let y3 = F::sub(&F::mul(&slope, &F::sub(left_x, &x3)), left_y);

        Ok(AffinePoint::new(x3, y3))
    }

    /// Doubles one already-validated finite affine point that is not
    /// `2`-torsion.
    ///
    /// The tangent slope is
    ///
    /// `λ = (3x^2 + 2Ax + 1) / (2By)`
    ///
    /// and the doubled point is again recovered from
    ///
    /// `x3 = B λ^2 - A - 2x`, `y3 = λ (x - x3) - y`.
    fn double_non_two_torsion_finite_point(
        &self,
        x: &F::Elem,
        y: &F::Elem,
    ) -> Result<AffinePoint<F>, CurveError> {
        let slope = F::div(
            &self.doubling_slope_numerator(x),
            &F::mul(&F::from_i64(2), &F::mul(self.b(), y)),
        )
        .expect("non-2-torsion affine points have nonzero tangent denominator");
        let x3 = F::sub(
            &F::sub(&F::mul(self.b(), &F::square(&slope)), self.a()),
            &F::mul(&F::from_i64(2), x),
        );
        let y3 = F::sub(&F::mul(&slope, &F::sub(x, &x3)), y);

        Ok(AffinePoint::new(x3, y3))
    }

    fn add_unchecked(
        &self,
        left: &AffinePoint<F>,
        right: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        match (left, right) {
            (AffinePoint::Infinity, _) => Ok(right.clone()),
            (_, AffinePoint::Infinity) => Ok(left.clone()),
            (
                AffinePoint::Finite {
                    x: left_x,
                    y: left_y,
                },
                AffinePoint::Finite {
                    x: right_x,
                    y: right_y,
                },
            ) => {
                if F::eq(left_x, right_x) {
                    if F::eq(left_y, &self.affine_negated_y(right_y)) {
                        return Ok(AffinePoint::Infinity);
                    }
                    self.double_non_two_torsion_finite_point(left_x, left_y)
                } else {
                    self.add_distinct_finite_points(left_x, left_y, right_x, right_y)
                }
            }
        }
    }
}

impl<F: Field> GroupCurveModel for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    fn neg(&self, point: &Self::Point) -> Self::Point {
        match point {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => {
                let negated = AffinePoint::new(x.clone(), self.affine_negated_y(y));
                debug_assert!(
                    !self.contains_affine_point(point) || self.contains_affine_point(&negated),
                    "Montgomery negation should preserve on-curve inputs"
                );
                negated
            }
        }
    }

    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }
        self.add_unchecked(left, right)
    }

    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }
        let negated = self.neg(right);
        self.add_unchecked(left, &negated)
    }

    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        self.add_unchecked(point, point)
    }

    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_unchecked(&result, &base)?;
            }
            k >>= 1;
            if k > 0 {
                base = self.double(&base)?;
            }
        }
        Ok(result)
    }

    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            Err(CurveError::PointNotOnCurve)
        } else if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar(point, scalar as u64)
        }
    }
}
