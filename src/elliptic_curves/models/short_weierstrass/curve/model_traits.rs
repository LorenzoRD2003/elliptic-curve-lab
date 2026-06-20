use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    traits::{
        AffineCurveModel, CurveModel, HasJInvariant, LiftXCoordinate, LiftedPoints,
        RelativeFrobeniusCurveModel,
    },
};
use crate::fields::traits::{Field, FiniteField, SqrtField};

impl<F: Field> CurveModel for ShortWeierstrassCurve<F> {
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
        match point {
            Self::Point::Infinity => true,
            Self::Point::Finite { x, y } => {
                let left = F::square(y);
                let right = self.rhs_value(x);
                F::eq(&left, &right)
            }
        }
    }
}

impl<F: Field> HasJInvariant for ShortWeierstrassCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        ShortWeierstrassCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for ShortWeierstrassCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = self.unchecked_point(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for ShortWeierstrassCurve<F> {
    fn lift_x(&self, x: Self::Elem) -> Result<LiftedPoints<Self::Point>, CurveError> {
        let (left_y, right_y) = match F::sqrt_pair(&self.rhs_value(&x)) {
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

impl<F: FiniteField> RelativeFrobeniusCurveModel for ShortWeierstrassCurve<F> {
    fn relative_frobenius(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        self.relative_frobenius_point(point)
    }
}
