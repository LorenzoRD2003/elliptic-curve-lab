use crate::elliptic_curves::{
    AffinePoint, CurveError, ProjectivePoint, ShortWeierstrassCurve,
    traits::{CurveModel, HasProjectiveModel},
};
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassCurve<F> {
    pub(super) fn contains_projective_point(&self, point: &ProjectivePoint<F>) -> bool {
        match point {
            ProjectivePoint::Infinity => true,
            ProjectivePoint::Finite { x, y, z } => {
                if F::is_zero(z) {
                    return false;
                }

                let left = F::mul(&F::square(y), z);
                let z2 = F::square(z);
                let right = F::add(
                    &F::add(&F::cube(x), &F::mul(self.a(), &F::mul(x, &z2))),
                    &F::mul(self.b(), &F::mul(&z2, z)),
                );
                F::eq(&left, &right)
            }
        }
    }

    pub(super) fn projective_point_from_affine(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ProjectivePoint<F>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        Ok(ProjectivePoint::from_affine(point))
    }

    pub(super) fn projective_point_to_affine(
        &self,
        point: &ProjectivePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        let affine = point.to_affine()?;
        if self.contains(&affine) {
            Ok(affine)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: Field> HasProjectiveModel for ShortWeierstrassCurve<F> {
    type ProjectivePoint = ProjectivePoint<F>;

    fn to_projective(&self, point: &Self::Point) -> Result<Self::ProjectivePoint, CurveError> {
        self.projective_point_from_affine(point)
    }

    fn to_affine_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::Point, CurveError> {
        self.projective_point_to_affine(point)
    }

    fn is_projective_point_on_curve(&self, point: &Self::ProjectivePoint) -> bool {
        self.contains_projective_point(point)
    }

    fn projective_identity(&self) -> Self::ProjectivePoint {
        ProjectivePoint::Infinity
    }
}
