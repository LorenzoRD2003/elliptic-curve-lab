use num_bigint::BigUint;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{HasseInterval, group_order::MestreSide},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(super) fn unique_mestre_multiple(
        &self,
        interval: &HasseInterval,
        lower_bound: &BigUint,
    ) -> Option<u128> {
        lower_bound
            .try_into()
            .ok()
            .and_then(|bound: u128| interval.unique_multiple_of(bound))
    }

    pub(super) fn curve_for_side<'a>(
        &'a self,
        side: MestreSide,
        twist_curve: &'a ShortWeierstrassCurve<F>,
    ) -> &'a ShortWeierstrassCurve<F> {
        match side {
            MestreSide::Original => self,
            MestreSide::QuadraticTwist => twist_curve,
        }
    }

    pub(super) fn lower_bound_for_side_mut<'a>(
        &self,
        side: MestreSide,
        original_lower_bound: &'a mut BigUint,
        twist_lower_bound: &'a mut BigUint,
    ) -> &'a mut BigUint {
        match side {
            MestreSide::Original => original_lower_bound,
            MestreSide::QuadraticTwist => twist_lower_bound,
        }
    }

    pub(super) fn unique_multiple_for_side(
        &self,
        interval: &HasseInterval,
        side: MestreSide,
        original_lower_bound: &BigUint,
        twist_lower_bound: &BigUint,
    ) -> Option<u128> {
        let lower_bound = match side {
            MestreSide::Original => original_lower_bound,
            MestreSide::QuadraticTwist => twist_lower_bound,
        };
        self.unique_mestre_multiple(interval, lower_bound)
    }

    pub(super) fn next_mestre_side(&self, side: MestreSide) -> MestreSide {
        match side {
            MestreSide::Original => MestreSide::QuadraticTwist,
            MestreSide::QuadraticTwist => MestreSide::Original,
        }
    }
}
