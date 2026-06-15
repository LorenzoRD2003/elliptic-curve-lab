use num_bigint::BigUint;

use super::{
    ExponentAccumulationReport, ExponentAccumulationStep, GroupExponentReport,
    GroupExponentStrategy,
};
use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    traits::{EnumerableCurveModel, FiniteGroupCurveModel, PointIndexSampler},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};
use crate::numerics::integer_arithmetic::lcm_biguint;

impl<F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    /// Recovers or estimates `λ(E(F_q))` by one requested strategy.
    ///
    /// Complexity:
    /// - [`GroupExponentStrategy::Exhaustive`]: `Θ(#E(F_q)^2)` group additions
    ///   in the current direct-traversal implementation, because it computes
    ///   every point order and then takes their `lcm`.
    /// - [`GroupExponentStrategy::RandomPoints`]: `Θ(s)` calls to
    ///   [`Self::point_order_by`] for `s = max_samples`, plus `Θ(s)` exact
    ///   `lcm` updates on arbitrary-precision integers.
    ///
    /// The random-point route samples with replacement because each call to
    /// [`crate::elliptic_curves::EnumerableCurveModel::random_point`] draws
    /// from a freshly materialized point list.
    pub fn group_exponent_by<S: PointIndexSampler>(
        &self,
        strategy: GroupExponentStrategy,
        sampler: &mut S,
    ) -> Result<GroupExponentReport<AffinePoint<F>>, CurveError> {
        match strategy {
            GroupExponentStrategy::Exhaustive => Ok(GroupExponentReport::Exhaustive(
                BigUint::from(self.exponent()),
            )),
            GroupExponentStrategy::RandomPoints {
                max_samples,
                point_order_strategy,
            } => {
                let mut steps = Vec::with_capacity(max_samples);
                let mut accumulated_lcm = BigUint::from(1u8);

                for _ in 0..max_samples {
                    let Some(point) = self.random_point(sampler) else {
                        break;
                    };
                    let point_order_report =
                        self.point_order_by(&point, point_order_strategy.clone())?;
                    accumulated_lcm =
                        lcm_biguint(&accumulated_lcm, point_order_report.exact_order());
                    steps.push(ExponentAccumulationStep::new(
                        point,
                        point_order_report,
                        accumulated_lcm.clone(),
                    ));
                }

                Ok(GroupExponentReport::RandomPoints(Box::new(
                    ExponentAccumulationReport::from_steps(
                        max_samples,
                        point_order_strategy,
                        steps,
                    ),
                )))
            }
        }
    }
}
