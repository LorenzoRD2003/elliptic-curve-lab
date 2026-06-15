use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusTrace, HasseInterval,
        group_order::GroupOrderReport,
        group_order::{MestreConfig, MestreGroupOrderReport, MestreSide, MestreStepReport},
    },
};
use crate::fields::{
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField},
};

pub(super) struct MestreFinalizeInput {
    pub(super) config: MestreConfig,
    pub(super) base_field: FiniteFieldDescriptor,
    pub(super) prime: u128,
    pub(super) interval: HasseInterval,
    pub(super) original_lower_bound: BigUint,
    pub(super) twist_lower_bound: BigUint,
    pub(super) steps: Vec<MestreStepReport>,
}

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(super) fn finalize_mestre_group_order_report(
        &self,
        input: MestreFinalizeInput,
    ) -> Result<GroupOrderReport, CurveError> {
        let (resolved_side, original_order, twist_order) = self.resolve_mestre_group_orders(
            &input.interval,
            input.prime,
            &input.original_lower_bound,
            &input.twist_lower_bound,
        );

        let original_trace = FrobeniusTrace::from_order(
            input.base_field.clone(),
            u64::try_from(original_order).expect("Mestre over u64 prime fields should fit in u64"),
        )?;
        let twist_trace = FrobeniusTrace::from_order(
            input.base_field,
            u64::try_from(twist_order).expect("Mestre over u64 prime fields should fit in u64"),
        )?;

        Ok(GroupOrderReport::MestreFp(Box::new(
            MestreGroupOrderReport::new(
                input.config,
                resolved_side,
                original_trace,
                twist_trace,
                input.steps,
            ),
        )))
    }

    fn resolve_mestre_group_orders(
        &self,
        interval: &HasseInterval,
        prime: u128,
        original_lower_bound: &BigUint,
        twist_lower_bound: &BigUint,
    ) -> (MestreSide, u128, u128) {
        if let Some(original_order) = self.unique_multiple_for_side(
            interval,
            MestreSide::Original,
            original_lower_bound,
            twist_lower_bound,
        ) {
            (
                MestreSide::Original,
                original_order,
                self.complementary_group_order_from_twist_sum(prime, original_order),
            )
        } else {
            let twist_order = self
                .unique_multiple_for_side(
                    interval,
                    MestreSide::QuadraticTwist,
                    original_lower_bound,
                    twist_lower_bound,
                )
                .expect("the loop exits only when one Mestre side is unique");
            (
                MestreSide::QuadraticTwist,
                self.complementary_group_order_from_twist_sum(prime, twist_order),
                twist_order,
            )
        }
    }

    fn complementary_group_order_from_twist_sum(&self, prime: u128, known_order: u128) -> u128 {
        prime
            .checked_mul(2)
            .and_then(|double_prime| double_prime.checked_add(2))
            .and_then(|order_sum| order_sum.checked_sub(known_order))
            .expect("2p + 2 - known twist-related order should stay in range")
    }
}
