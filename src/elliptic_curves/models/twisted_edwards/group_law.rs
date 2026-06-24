use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    traits::{CurveModel, GroupCurveModel},
};
use crate::fields::{FieldError, traits::Field};

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn affine_neg_x(&self, x: &F::Elem) -> F::Elem {
        F::neg(x)
    }

    fn addition_x_numerator(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> F::Elem {
        F::add(&F::mul(left_x, right_y), &F::mul(left_y, right_x))
    }

    fn addition_x_denominator(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> F::Elem {
        F::add(
            &F::one(),
            &F::mul(
                self.d(),
                &F::mul(&F::mul(left_x, right_x), &F::mul(left_y, right_y)),
            ),
        )
    }

    fn addition_y_numerator(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> F::Elem {
        F::sub(
            &F::mul(left_y, right_y),
            &F::mul(self.a(), &F::mul(left_x, right_x)),
        )
    }

    fn addition_y_denominator(
        &self,
        left_x: &F::Elem,
        left_y: &F::Elem,
        right_x: &F::Elem,
        right_y: &F::Elem,
    ) -> F::Elem {
        F::sub(
            &F::one(),
            &F::mul(
                self.d(),
                &F::mul(&F::mul(left_x, right_x), &F::mul(left_y, right_y)),
            ),
        )
    }

    /// Adds two already-validated affine points with the generic twisted-Edwards
    /// formulas
    ///
    /// `x3 = (x1 y2 + y1 x2) / (1 + d x1 x2 y1 y2)`
    ///
    /// `y3 = (y1 y2 - a x1 x2) / (1 - d x1 x2 y1 y2)`.
    ///
    /// This staged affine implementation is honest about the fact that these
    /// denominators can vanish on valid inputs for generic `(a, d)`. In that
    /// case it returns `CurveError::Field(FieldError::DivisionByZero)`.
    fn add_unchecked(
        &self,
        left: &AffinePoint<F>,
        right: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if self.is_identity(left) {
            return Ok(right.clone());
        }
        if self.is_identity(right) {
            return Ok(left.clone());
        }

        let (
            AffinePoint::Finite {
                x: left_x,
                y: left_y,
            },
            AffinePoint::Finite {
                x: right_x,
                y: right_y,
            },
        ) = (left, right)
        else {
            unreachable!("non-identity twisted-Edwards points are always finite");
        };

        if F::eq(left_x, &self.affine_neg_x(right_x)) && F::eq(left_y, right_y) {
            return Ok(self.identity());
        }

        let x_denominator = self.addition_x_denominator(left_x, left_y, right_x, right_y);
        let y_denominator = self.addition_y_denominator(left_x, left_y, right_x, right_y);

        if F::is_zero(&x_denominator) || F::is_zero(&y_denominator) {
            return Err(CurveError::Field(FieldError::DivisionByZero));
        }

        let x3 = F::div(
            &self.addition_x_numerator(left_x, left_y, right_x, right_y),
            &x_denominator,
        )
        .expect("checked non-zero twisted-Edwards x denominator");
        let y3 = F::div(
            &self.addition_y_numerator(left_x, left_y, right_x, right_y),
            &y_denominator,
        )
        .expect("checked non-zero twisted-Edwards y denominator");

        Ok(AffinePoint::new(x3, y3))
    }

    fn double_unchecked(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, CurveError> {
        self.add_unchecked(point, point)
    }

    fn mul_scalar_unchecked(
        &self,
        point: &AffinePoint<F>,
        scalar: u64,
    ) -> Result<AffinePoint<F>, CurveError> {
        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_unchecked(&result, &base)?;
            }
            k >>= 1;
            if k > 0 {
                base = self.double_unchecked(&base)?;
            }
        }

        Ok(result)
    }
}

impl<F: Field> GroupCurveModel for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn neg(&self, point: &Self::Point) -> Self::Point {
        match point {
            AffinePoint::Infinity => AffinePoint::Infinity,
            AffinePoint::Finite { x, y } => {
                let negated = AffinePoint::new(self.affine_neg_x(x), y.clone());
                debug_assert!(
                    !self.contains_affine_point(point) || self.contains_affine_point(&negated),
                    "twisted-Edwards negation should preserve on-curve inputs"
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

        self.double_unchecked(point)
    }

    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.mul_scalar_unchecked(point, scalar)
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
            self.mul_scalar_unchecked(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar_unchecked(point, scalar as u64)
        }
    }
}
