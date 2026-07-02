/// Curve models equipped with an explicit additive group law on their points.
///
/// This trait is narrower than [`CurveModel`]: it is for concrete models where
/// the crate is ready to expose actual point addition, doubling, and scalar
/// multiplication, not just curve membership and constructors.
///
/// Implementations should treat `CurveError::PointNotOnCurve` as the honest
/// failure mode for finite inputs that do not belong to the model.
pub trait GroupCurveModel: CurveModel
where
    Self::Point: Clone,
{
    /// Returns the additive inverse of a point.
    ///
    /// For on-curve inputs this should stay on the curve and satisfy
    /// `P + (-P) = O`.
    fn neg(&self, point: &Self::Point) -> Self::Point;

    /// Adds two curve points.
    ///
    /// Implementations should return the honest group-law sum when both inputs
    /// lie on the curve, including the usual identity and inverse cases.
    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError>;

    /// Subtracts `right` from `left`.
    ///
    /// The default implementation negates `right` and then reuses
    /// [`Self::add`].
    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add(left, &negated)
    }

    /// Doubles a point under the group law.
    ///
    /// The default implementation calls [`Self::add`] with the same point
    /// twice after validating membership.
    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.add(point, point)
    }

    /// Multiplies a point by a non-negative integer using double-and-add.
    ///
    /// This is the clear educational baseline rather than an optimized
    /// cryptographic ladder.
    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add(&result, &base)?;
            }

            k >>= 1;

            if k > 0 {
                base = self.double(&base)?;
            }
        }

        Ok(result)
    }

    /// Multiplies a point by a signed integer.
    ///
    /// Negative scalars are handled by negating the point and reusing
    /// [`Self::mul_scalar`].
    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar(point, scalar as u64)
        }
    }

    /// Returns whether the point is killed by multiplication by `n`.
    ///
    /// This helper is meant for the usual positive-integer notion of
    /// `n`-torsion. To avoid the degenerate convention `[0]P = O`, this method
    /// returns `false` when `n == 0`.
    ///
    /// Invalid off-curve inputs are treated honestly and return `false`.
    fn is_torsion_point(&self, point: &Self::Point, n: u64) -> bool {
        if n == 0 {
            return false;
        }

        self.mul_scalar(point, n)
            .is_ok_and(|multiple| self.is_identity(&multiple))
    }
}
use crate::elliptic_curves::{CurveError, traits::CurveModel};
