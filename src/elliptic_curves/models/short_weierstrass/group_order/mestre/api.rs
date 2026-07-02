use std::hash::Hash;

use super::{finalize::MestreFinalizeInput, loop_logic::MestreLoopState};
use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{HasseInterval, group_order::GroupOrderReport, group_order::MestreConfig},
    traits::PointIndexSampler,
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(crate) fn group_order_by_mestre_fp<S: PointIndexSampler>(
        &self,
        config: MestreConfig,
        sampler: &mut S,
    ) -> Result<GroupOrderReport, CurveError>
    where
        F::Elem: Hash,
    {
        let (base_field, prime) = self.validate_mestre_prime_field()?;
        let interval = HasseInterval::for_q(&prime)?;
        let twist_curve = self.select_genuine_quadratic_twist_for_mestre()?;
        let mut state = MestreLoopState::new();

        while self.mestre_needs_more_steps(&interval, &state) {
            if config
                .max_iterations()
                .is_some_and(|max_iterations| state.iterations >= max_iterations)
            {
                return Err(CurveError::MestreIterationCapReached {
                    max_iterations: config
                        .max_iterations()
                        .expect("checked presence above for the capped branch"),
                });
            }

            self.run_one_mestre_step(&interval, &twist_curve, &mut state, sampler)?;
        }

        self.finalize_mestre_group_order_report(MestreFinalizeInput {
            config,
            base_field,
            prime,
            interval,
            original_lower_bound: state.original_lower_bound,
            twist_lower_bound: state.twist_lower_bound,
            steps: state.steps,
        })
    }
}
