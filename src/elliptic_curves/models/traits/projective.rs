use core::fmt;

use crate::elliptic_curves::{
    CurveError,
    traits::{CurveModel, ScalarInput},
};
use num_bigint::BigUint;
use num_traits::{One, Zero};

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
    ///
    /// The default implementation is a correctness baseline: it lifts the
    /// affine input into projective form and delegates to `add_projective`.
    /// Models may override this with a genuinely specialized mixed-add formula,
    /// for example one that exploits `Z_2 = 1`.
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
        scalar: impl ScalarInput,
    ) -> Result<Self::ProjectivePoint, CurveError> {
        if !self.is_projective_point_on_curve(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut k = scalar.into_biguint_scalar();
        let mut result = self.projective_identity();
        let mut base = point.clone();

        while !k.is_zero() {
            if (&k & BigUint::one()) == BigUint::one() {
                result = self.add_projective(&result, &base)?;
            }
            k >>= 1usize;
            if !k.is_zero() {
                base = self.double_projective(&base)?;
            }
        }
        Ok(result)
    }
}
