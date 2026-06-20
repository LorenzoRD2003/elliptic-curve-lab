use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve,
    traits::{CurveModel, GroupCurveModel},
};
use crate::fields::traits::Field;

impl<F: Field> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    // TODO: Promote this native affine implementation to projective coordinates.
    // The current formulas are mathematically honest and remove the old
    // dependency on short-Weierstrass reduction, but they still branch on
    // affine denominators and therefore are not the long-term execution model
    // we want for the general family.

    fn negated_y_coordinate(&self, x: &F::Elem, y: &F::Elem) -> F::Elem {
        F::sub(&F::neg(y), &F::add(&F::mul(self.a1(), x), self.a3()))
    }

    fn affine_sum_from_line(
        &self,
        x1: &F::Elem,
        x2: &F::Elem,
        y1: &F::Elem,
        slope: &F::Elem,
    ) -> AffinePoint<F> {
        let intercept = F::sub(y1, &F::mul(slope, x1));
        let x3 = F::sub(
            &F::sub(
                &F::add(&F::square(slope), &F::mul(self.a1(), slope)),
                self.a2(),
            ),
            &F::add(x1, x2),
        );
        let y3 = F::sub(
            &F::neg(&F::add(&F::mul(&F::add(slope, self.a1()), &x3), &intercept)),
            self.a3(),
        );
        let sum = AffinePoint::new(x3, y3);
        debug_assert!(self.contains_affine_point(&sum));
        sum
    }

    fn affine_add_finite_points(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> Result<AffinePoint<F>, CurveError> {
        if F::eq(left_x, right_x) {
            if F::eq(left_y, right_y) {
                return self.affine_double_finite_point(left_x, left_y);
            }

            return Ok(self.identity());
        }

        let slope = F::div(&F::sub(right_y, left_y), &F::sub(right_x, left_x))?;
        Ok(self.affine_sum_from_line(left_x, right_x, left_y, &slope))
    }

    fn affine_double_finite_point(
        &self,
        x: &F::Elem,
        y: &F::Elem,
    ) -> Result<AffinePoint<F>, CurveError> {
        let denominator = F::add(
            &F::mul(&F::from_i64(2), y),
            &F::add(&F::mul(self.a1(), x), self.a3()),
        );

        if F::is_zero(&denominator) {
            return Ok(self.identity());
        }

        let numerator = F::sub(
            &F::add(
                &F::add(
                    &F::mul(&F::from_i64(3), &F::square(x)),
                    &F::mul(&F::mul(&F::from_i64(2), self.a2()), x),
                ),
                self.a4(),
            ),
            &F::mul(self.a1(), y),
        );

        let slope = F::div(&numerator, &denominator)?;
        Ok(self.affine_sum_from_line(x, x, y, &slope))
    }
}

impl<F: Field> GroupCurveModel for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn neg(&self, point: &Self::Point) -> Self::Point {
        match point {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => {
                let negated = AffinePoint::new(x.clone(), self.negated_y_coordinate(x, y));
                debug_assert!(
                    !self.contains_affine_point(point) || self.contains_affine_point(&negated),
                    "general-Weierstrass negation should preserve on-curve inputs"
                );
                negated
            }
        }
    }

    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        match (left, right) {
            (AffinePoint::Infinity, point) | (point, AffinePoint::Infinity) => Ok(point.clone()),
            (
                AffinePoint::Finite {
                    x: left_x,
                    y: left_y,
                },
                AffinePoint::Finite {
                    x: right_x,
                    y: right_y,
                },
            ) => self.affine_add_finite_points(left_x, left_y, right_x, right_y),
        }
    }

    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(left) || !self.contains_affine_point(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add(left, &negated)
    }

    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => self.affine_double_finite_point(x, y),
        }
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
                result = self.add(&result, &base)?;
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
            return Err(CurveError::PointNotOnCurve);
        }

        if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar(point, scalar as u64)
        }
    }
}
