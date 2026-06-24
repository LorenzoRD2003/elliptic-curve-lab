use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints},
};
use crate::fields::traits::{Field, SqrtField};

impl<F: Field> CurveModel for TwistedEdwardsCurve<F> {
    type Elem = F::Elem;
    type BaseField = F;
    type Point = AffinePoint<F>;

    fn identity(&self) -> Self::Point {
        AffinePoint::new(F::zero(), F::one())
    }

    fn is_identity(&self, point: &Self::Point) -> bool {
        match point {
            AffinePoint::Infinity => false,
            AffinePoint::Finite { x, y } => F::is_zero(x) && F::eq(y, &F::one()),
        }
    }

    fn contains(&self, point: &Self::Point) -> bool {
        self.contains_affine_point(point)
    }
}

impl<F: Field> HasJInvariant for TwistedEdwardsCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        TwistedEdwardsCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for TwistedEdwardsCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = AffinePoint::new(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for TwistedEdwardsCurve<F> {
    /// Recovers the affine fiber above one chosen `x` on the twisted-Edwards
    /// model `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2`.
    ///
    /// Rearranging for one fixed `x` gives the lifting problem
    ///
    /// `y^2 = (1 - a x^2) / (1 - d x^2)`
    ///
    /// whenever `1 - d x^2 ≠ 0`. Because validated twisted-Edwards curves
    /// satisfy `a ≠ d`, the numerator and denominator cannot vanish
    /// simultaneously, so a vanishing denominator means the fiber is empty.
    fn lift_x(&self, x: Self::Elem) -> Result<LiftedPoints<Self::Point>, CurveError> {
        let denominator = self.y_sq_denominator(&x);
        if F::is_zero(&denominator) {
            return Ok(LiftedPoints::NoPoint);
        }

        let y_sq = F::div(&self.y_sq_numerator(&x), &denominator)
            .expect("checked non-zero twisted-Edwards lifting denominator");
        let (left_y, right_y) = match F::sqrt_pair(&y_sq) {
            Some(pair) => pair,
            None => return Ok(LiftedPoints::NoPoint),
        };

        let left = self.point(x.clone(), left_y)?;
        let right = self.point(x, right_y)?;

        if left == right {
            Ok(LiftedPoints::OnePoint(left))
        } else {
            Ok(LiftedPoints::TwoPoints(left, right))
        }
    }
}

impl<F: Field> TwistedEdwardsCurve<F> {
    fn y_sq_numerator(&self, x: &F::Elem) -> F::Elem {
        F::sub(&F::one(), &F::mul(self.a(), &F::square(x)))
    }

    fn y_sq_denominator(&self, x: &F::Elem) -> F::Elem {
        F::sub(&F::one(), &F::mul(self.d(), &F::square(x)))
    }
}
