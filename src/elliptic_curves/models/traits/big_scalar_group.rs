use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::elliptic_curves::{CurveError, traits::GroupCurveModel};

/// Curve models that support internal `BigUint` scalar multiplication.
///
/// This trait refines [`GroupCurveModel`] with a default educational helper
/// used by algorithms that already carry arbitrary-precision scalar data.
pub trait BigScalarGroupCurveModel: GroupCurveModel {
    /// Multiplies one curve point by a non-negative `BigUint` scalar.
    ///
    /// This is the internal big-integer analogue of the public scalar-input
    /// double-and-add surface on [`GroupCurveModel`].
    fn mul_scalar_biguint(
        &self,
        point: &Self::Point,
        scalar: &BigUint,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar.clone();

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
}

impl<T: GroupCurveModel + ?Sized> BigScalarGroupCurveModel for T {}
