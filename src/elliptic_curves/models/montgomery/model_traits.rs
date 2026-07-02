use crate::elliptic_curves::{
    AffinePoint, CurveError, MontgomeryCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints},
};
use crate::fields::traits::SqrtField;
use crate::fields::traits::*;

impl<F: Field> CurveModel for MontgomeryCurve<F> {
    type Elem = F::Elem;
    type BaseField = F;
    type Point = AffinePoint<F>;

    fn identity(&self) -> Self::Point {
        Self::Point::infinity()
    }

    fn is_identity(&self, point: &Self::Point) -> bool {
        point.is_identity()
    }

    fn contains(&self, point: &Self::Point) -> bool {
        self.contains_affine_point(point)
    }
}

impl<F: Field> HasJInvariant for MontgomeryCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        MontgomeryCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for MontgomeryCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = AffinePoint::new(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for MontgomeryCurve<F> {
    /// Recovers the affine fiber above one chosen `x` on the Montgomery model
    ///
    /// `B y^2 = x^3 + A x^2 + x`.
    ///
    /// Since validated Montgomery curves satisfy `B != 0`, the lifting problem
    /// is the square-root equation
    ///
    /// `y^2 = (x^3 + A x^2 + x) / B`.
    fn lift_x(&self, x: Self::Elem) -> Result<LiftedPoints<Self::Point>, CurveError> {
        let y_squared =
            F::div(&self.rhs_value(&x), self.b()).expect("validated Montgomery curve has B != 0");
        let (left_y, right_y) = match F::sqrt_pair(&y_squared) {
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
