use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve, traits::CurveModel};
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
