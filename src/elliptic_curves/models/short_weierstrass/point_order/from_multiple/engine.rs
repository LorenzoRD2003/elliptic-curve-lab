use num_bigint::BigUint;
use num_traits::One;

use super::{PointOrderFromMultipleReport, PointOrderReductionStep};
use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    group_algorithms::CyclicPrimaryOrderGroupCurveModel, traits::BigScalarGroupCurveModel,
};
use crate::fields::traits::FiniteField;
use crate::numerics::PrimePowerTable;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    pub(super) fn recover_point_order_from_normalized_factorization(
        &self,
        point: &AffinePoint<F>,
        supplied_multiple: BigUint,
        normalized_factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        let mut remaining_multiple = supplied_multiple.clone();
        let mut exact_order = BigUint::one();
        let mut steps = Vec::with_capacity(normalized_factorization.len());

        for (prime, exponent_in_multiple) in normalized_factorization {
            let powers = PrimePowerTable::up_through(prime, *exponent_in_multiple);
            let prime_power = powers.power(*exponent_in_multiple);
            let cofactor = &remaining_multiple / prime_power;
            let primary_component = if cofactor == BigUint::one() {
                point.clone()
            } else {
                self.mul_scalar_biguint(point, &cofactor)?
            };

            let local_report = self.recover_cyclic_primary_order(&primary_component, &powers)?;
            let removed_exponent = local_report.removed_exponent();
            let local_exact_power = powers.power(local_report.exact_exponent());
            let removed_power = powers.power(removed_exponent);

            exact_order *= local_exact_power;
            remaining_multiple /= removed_power;

            steps.push(PointOrderReductionStep::new(
                local_report.prime().clone(),
                local_report.exponent_bound(),
                removed_exponent,
                remaining_multiple.clone(),
            ));
        }

        Ok(PointOrderFromMultipleReport::new(
            supplied_multiple,
            exact_order,
            steps,
        ))
    }
}
