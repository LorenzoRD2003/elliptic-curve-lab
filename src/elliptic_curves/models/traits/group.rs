use crate::elliptic_curves::{
    CurveError,
    traits::{CurveModel, ScalarInput},
};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};

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
    fn mul_scalar(
        &self,
        point: &Self::Point,
        scalar: impl ScalarInput,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut k = scalar.into_biguint_scalar();
        let mut result = self.identity();
        let mut base = point.clone();

        while !k.is_zero() {
            if (&k & BigUint::one()) == BigUint::one() {
                result = self.add(&result, &base)?;
            }

            k >>= 1usize;

            if !k.is_zero() {
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
        scalar: impl Into<BigInt>,
    ) -> Result<Self::Point, CurveError> {
        let scalar = scalar.into();
        if scalar.sign() == Sign::Minus {
            let negated = self.neg(point);
            self.mul_scalar(
                &negated,
                (-scalar)
                    .to_biguint()
                    .expect("negated negative scalar should be non-negative"),
            )
        } else {
            self.mul_scalar(
                point,
                scalar
                    .to_biguint()
                    .expect("non-negative scalar should convert to BigUint"),
            )
        }
    }

    /// Returns whether the point is killed by multiplication by `n`.
    ///
    /// This helper is meant for the usual positive-integer notion of
    /// `n`-torsion. To avoid the degenerate convention `[0]P = O`, this method
    /// returns `false` when `n == 0`.
    ///
    /// Invalid off-curve inputs are treated honestly and return `false`.
    fn is_torsion_point(&self, point: &Self::Point, n: impl ScalarInput) -> bool {
        let n = n.into_biguint_scalar();
        if n.is_zero() {
            return false;
        }

        self.mul_scalar(point, n)
            .is_ok_and(|multiple| self.is_identity(&multiple))
    }
}
