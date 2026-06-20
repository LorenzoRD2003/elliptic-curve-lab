use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant},
};
use crate::fields::traits::Field;

impl<F: Field> CurveModel for GeneralWeierstrassCurve<F> {
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

impl<F: Field> HasJInvariant for GeneralWeierstrassCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        GeneralWeierstrassCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for GeneralWeierstrassCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = AffinePoint::new(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}
