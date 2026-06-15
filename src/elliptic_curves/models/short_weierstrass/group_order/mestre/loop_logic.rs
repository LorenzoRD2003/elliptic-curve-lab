use num_bigint::BigUint;
use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{HasseInterval, group_order::MestreSide, group_order::MestreStepReport},
    traits::{EnumerableCurveModel, PointIndexSampler},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};
use crate::numerics::{NormalizedPrimePowerFactorization, integer_arithmetic::lcm_biguint};

pub(super) struct MestreLoopState {
    pub(super) original_lower_bound: BigUint,
    pub(super) twist_lower_bound: BigUint,
    pub(super) next_side: MestreSide,
    pub(super) steps: Vec<MestreStepReport>,
    pub(super) iterations: usize,
}

impl MestreLoopState {
    pub(super) fn new() -> Self {
        Self {
            original_lower_bound: BigUint::from(1u8),
            twist_lower_bound: BigUint::from(1u8),
            next_side: MestreSide::Original,
            steps: Vec::new(),
            iterations: 0,
        }
    }
}

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(super) fn mestre_side_is_unique(
        &self,
        interval: &HasseInterval,
        side: MestreSide,
        state: &MestreLoopState,
    ) -> bool {
        self.unique_multiple_for_side(
            interval,
            side,
            &state.original_lower_bound,
            &state.twist_lower_bound,
        )
        .is_some()
    }

    pub(super) fn mestre_needs_more_steps(
        &self,
        interval: &HasseInterval,
        state: &MestreLoopState,
    ) -> bool {
        !self.mestre_side_is_unique(interval, MestreSide::Original, state)
            && !self.mestre_side_is_unique(interval, MestreSide::QuadraticTwist, state)
    }

    pub(super) fn run_one_mestre_step<S: PointIndexSampler>(
        &self,
        interval: &HasseInterval,
        twist_curve: &ShortWeierstrassCurve<F>,
        state: &mut MestreLoopState,
        sampler: &mut S,
    ) -> Result<(), CurveError>
    where
        F::Elem: Hash,
    {
        let curve = self.curve_for_side(state.next_side, twist_curve);
        let lower_bound = self.lower_bound_for_side_mut(
            state.next_side,
            &mut state.original_lower_bound,
            &mut state.twist_lower_bound,
        );

        let point = curve
            .random_point(sampler)
            .ok_or(CurveError::MestreSamplerExhausted)?;
        let Some(annihilating_multiple) =
            curve.find_annihilating_multiple_in_interval_bsgs(&point, interval.clone())?
        else {
            return Err(CurveError::NoAnnihilatingMultipleInHasseInterval {
                lower: interval.lower(),
                upper: interval.upper(),
            });
        };

        let annihilating_multiple_biguint = BigUint::from(annihilating_multiple);
        let factorization =
            NormalizedPrimePowerFactorization::factor(&annihilating_multiple_biguint)
                .expect("an annihilating multiple in H(p) should admit a prime factorization")
                .into_factors();
        let point_order_report = curve.point_order_from_multiple_with_trusted_factorization(
            &point,
            annihilating_multiple_biguint,
            &factorization,
        )?;

        *lower_bound = lcm_biguint(lower_bound, point_order_report.exact_order());
        state.steps.push(MestreStepReport::new(
            state.next_side,
            annihilating_multiple,
            point_order_report,
            lower_bound.clone(),
        ));

        state.next_side = self.next_mestre_side(state.next_side);
        state.iterations += 1;
        Ok(())
    }
}
