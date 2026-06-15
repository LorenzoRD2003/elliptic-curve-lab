use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    traits::{BigScalarGroupCurveModel, CurveModel},
};
use crate::fields::traits::FiniteField;
use crate::numerics::NormalizedPrimePowerFactorization;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    pub(super) fn validate_point_order_from_multiple_inputs(
        &self,
        point: &AffinePoint<F>,
        multiple: &BigUint,
    ) -> Result<(), CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        if multiple.is_zero() {
            return Err(CurveError::InvalidPointOrderMultiple {
                multiple: multiple.clone(),
            });
        }
        self.ensure_multiple_annihilates_point(point, multiple)
    }

    pub(super) fn annihilates_multiple(
        &self,
        point: &AffinePoint<F>,
        multiple: &BigUint,
    ) -> Result<bool, CurveError> {
        let image = self.mul_scalar_biguint(point, multiple)?;
        Ok(self.is_identity(&image))
    }

    pub(super) fn ensure_multiple_annihilates_point(
        &self,
        point: &AffinePoint<F>,
        multiple: &BigUint,
    ) -> Result<(), CurveError> {
        if self.annihilates_multiple(point, multiple)? {
            Ok(())
        } else {
            Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
                multiple: multiple.clone(),
            })
        }
    }
}

pub(super) fn normalized_factorization(
    multiple: &BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<NormalizedPrimePowerFactorization, CurveError> {
    NormalizedPrimePowerFactorization::checked(multiple, factorization).map_err(|_| {
        CurveError::InvalidPointOrderMultipleFactorization {
            multiple: multiple.clone(),
        }
    })
}

pub(super) fn trusted_normalized_factorization(
    multiple: &BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<NormalizedPrimePowerFactorization, CurveError> {
    NormalizedPrimePowerFactorization::trusted(multiple, factorization).map_err(|_| {
        CurveError::InvalidPointOrderMultipleFactorization {
            multiple: multiple.clone(),
        }
    })
}
