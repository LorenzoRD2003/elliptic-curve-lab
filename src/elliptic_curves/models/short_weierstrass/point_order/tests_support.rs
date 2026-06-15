use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    models::short_weierstrass::point_order::{
        PointOrderFromMultipleReport, PointOrderReductionStep,
    },
    traits::{BigScalarGroupCurveModel, CurveModel},
};
use crate::fields::traits::FiniteField;
use crate::numerics::NormalizedPrimePowerFactorization;

pub(crate) fn point_order_from_multiple_baseline<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    multiple: BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<PointOrderFromMultipleReport, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }
    if multiple.is_zero() {
        return Err(CurveError::InvalidPointOrderMultiple {
            multiple: multiple.clone(),
        });
    }
    if !curve.is_identity(&curve.mul_scalar_biguint(point, &multiple)?) {
        return Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
            multiple: multiple.clone(),
        });
    }

    let normalized_factorization =
        NormalizedPrimePowerFactorization::checked(&multiple, factorization)
            .map_err(|_| CurveError::InvalidPointOrderMultipleFactorization {
                multiple: multiple.clone(),
            })?
            .into_factors();

    let supplied_multiple = multiple;
    let mut remaining_multiple = supplied_multiple.clone();
    let mut steps = Vec::with_capacity(normalized_factorization.len());

    for (prime, exponent_in_multiple) in &normalized_factorization {
        let mut removed_exponent = 0;

        for _ in 0..*exponent_in_multiple {
            let candidate_multiple = &remaining_multiple / prime;
            if curve.is_identity(&curve.mul_scalar_biguint(point, &candidate_multiple)?) {
                remaining_multiple = candidate_multiple;
                removed_exponent += 1;
            } else {
                break;
            }
        }

        steps.push(PointOrderReductionStep::new(
            prime.clone(),
            *exponent_in_multiple,
            removed_exponent,
            remaining_multiple.clone(),
        ));
    }

    Ok(PointOrderFromMultipleReport::new(
        supplied_multiple,
        remaining_multiple,
        steps,
    ))
}
