use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    traits::{CurveModel, HasProjectiveModel},
    twisted_edwards::projective::ExtendedTwistedEdwardsPoint,
};
use crate::fields::traits::Field;

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    pub(crate) fn extended_point_from_affine(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ExtendedTwistedEdwardsPoint<F>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        ExtendedTwistedEdwardsPoint::from_affine(point)
    }

    pub(crate) fn extended_point_to_affine(
        &self,
        point: &ExtendedTwistedEdwardsPoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.contains_extended_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        let affine = point.to_affine()?;
        if self.contains(&affine) {
            Ok(affine)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: Field> HasProjectiveModel for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    type ProjectivePoint = ExtendedTwistedEdwardsPoint<F>;

    fn to_projective(&self, point: &Self::Point) -> Result<Self::ProjectivePoint, CurveError> {
        self.extended_point_from_affine(point)
    }

    fn to_affine_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::Point, CurveError> {
        self.extended_point_to_affine(point)
    }

    fn is_projective_point_on_curve(&self, point: &Self::ProjectivePoint) -> bool {
        self.contains_extended_point(point)
    }

    fn projective_identity(&self) -> Self::ProjectivePoint {
        ExtendedTwistedEdwardsPoint::identity()
    }
}
