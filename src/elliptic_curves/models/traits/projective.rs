use core::fmt;

use crate::elliptic_curves::{CurveError, traits::CurveModel};

/// Curve models that admit an explicit projective representation of their
/// public point surface.
pub trait HasProjectiveModel: CurveModel {
    type ProjectivePoint: Clone + fmt::Debug + PartialEq;

    /// Lifts one public curve point into the model's projective representation.
    fn to_projective(&self, point: &Self::Point) -> Result<Self::ProjectivePoint, CurveError>;

    /// Recovers one public curve point from a projective representative.
    fn to_affine_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::Point, CurveError>;

    /// Returns whether a projective point belongs to the curve model.
    fn is_projective_point_on_curve(&self, point: &Self::ProjectivePoint) -> bool;

    /// Returns the distinguished identity element in projective form.
    fn projective_identity(&self) -> Self::ProjectivePoint;
}

/// Curve models with an explicit projective group-law surface.
pub trait ProjectiveGroupCurveModel: HasProjectiveModel {
    /// Returns the additive inverse of a projective point.
    fn neg_projective(&self, point: &Self::ProjectivePoint) -> Self::ProjectivePoint;

    /// Adds two projective points.
    fn add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError>;

    /// Doubles one projective point.
    fn double_projective(
        &self,
        point: &Self::ProjectivePoint,
    ) -> Result<Self::ProjectivePoint, CurveError>;

    /// Adds a projective point to one point in the model's canonical public
    /// representation.
    fn mixed_add_projective(
        &self,
        left: &Self::ProjectivePoint,
        right: &Self::Point,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        let projective_right = self.to_projective(right)?;
        self.add_projective(left, &projective_right)
    }

    /// Multiplies a projective point by a non-negative integer.
    fn mul_scalar_projective(
        &self,
        point: &Self::ProjectivePoint,
        scalar: u64,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.projective_identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add_projective(&result, &base)?;
            }
            k >>= 1;
            if k > 0 {
                base = self.double_projective(&base)?;
            }
        }
        Ok(result)
    }
}
