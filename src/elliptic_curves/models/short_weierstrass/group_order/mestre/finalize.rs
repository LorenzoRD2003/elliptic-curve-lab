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
    pub(super) prime: BigUint,
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
            &input.prime,
            &input.original_lower_bound,
            &input.twist_lower_bound,
        );

        let original_trace = FrobeniusTrace::from_order(input.base_field.clone(), original_order)?;
        let twist_trace = FrobeniusTrace::from_order(input.base_field, twist_order)?;

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
        prime: &BigUint,
        original_lower_bound: &BigUint,
        twist_lower_bound: &BigUint,
    ) -> (MestreSide, BigUint, BigUint) {
        if let Some(original_order) = self.unique_multiple_for_side(
            interval,
            MestreSide::Original,
            original_lower_bound,
            twist_lower_bound,
        ) {
            (
                MestreSide::Original,
                original_order.clone(),
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
                self.complementary_group_order_from_twist_sum(prime, twist_order.clone()),
                twist_order,
            )
        }
    }

    fn complementary_group_order_from_twist_sum(
        &self,
        prime: &BigUint,
        known_order: BigUint,
    ) -> BigUint {
        prime * BigUint::from(2u8) + BigUint::from(2u8) - known_order
    }
}
