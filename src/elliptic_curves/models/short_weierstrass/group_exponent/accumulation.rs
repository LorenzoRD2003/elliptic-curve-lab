use num_bigint::BigUint;

use crate::elliptic_curves::short_weierstrass::{
    group_exponent::GroupExponentStrategy,
    point_order::{PointOrderReport, PointOrderStrategy},
};

/// One random-point step in the running `lcm` accumulation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentAccumulationStep<P> {
    point: P,
    point_order_report: PointOrderReport<P>,
    accumulated_lcm: BigUint,
}

impl<P> ExponentAccumulationStep<P> {
    pub(crate) fn new(
        point: P,
        point_order_report: PointOrderReport<P>,
        accumulated_lcm: BigUint,
    ) -> Self {
        Self {
            point,
            point_order_report,
            accumulated_lcm,
        }
    }

    /// Returns the sampled point.
    pub fn point(&self) -> &P {
        &self.point
    }

    /// Returns the point-order report used for this sampled point.
    pub fn point_order_report(&self) -> &PointOrderReport<P> {
        &self.point_order_report
    }

    /// Returns the running `lcm` after processing this point.
    pub fn accumulated_lcm(&self) -> &BigUint {
        &self.accumulated_lcm
    }
}

/// Report for the random-point accumulation route to an exponent lower bound.
///
/// This route is heuristic: the final accumulated value is always a lower
/// bound for the true exponent, but it becomes exact only if the sampled point
/// orders already capture all prime-power factors of `λ(E(F_q))`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentAccumulationReport<P> {
    samples_requested: usize,
    point_order_strategy: PointOrderStrategy,
    exponent_lower_bound: BigUint,
    steps: Vec<ExponentAccumulationStep<P>>,
}

impl<P> ExponentAccumulationReport<P> {
    pub(crate) fn from_steps(
        samples_requested: usize,
        point_order_strategy: PointOrderStrategy,
        steps: Vec<ExponentAccumulationStep<P>>,
    ) -> Self {
        let exponent_lower_bound = steps
            .last()
            .map(|step| step.accumulated_lcm.clone())
            .unwrap_or_else(|| BigUint::from(1u8));

        Self {
            samples_requested,
            point_order_strategy,
            exponent_lower_bound,
            steps,
        }
    }

    /// Returns how many samples were requested.
    pub fn samples_requested(&self) -> usize {
        self.samples_requested
    }

    /// Returns how many samples were actually processed.
    pub fn samples_taken(&self) -> usize {
        self.steps.len()
    }

    /// Returns whether the sampler supplied all requested points.
    pub fn completed_requested_samples(&self) -> bool {
        self.samples_taken() == self.samples_requested
    }

    /// Returns the point-order strategy used on each sample.
    pub fn point_order_strategy(&self) -> &PointOrderStrategy {
        &self.point_order_strategy
    }

    /// Returns the accumulated lower bound for `λ(E(F_q))`.
    pub fn exponent_lower_bound(&self) -> &BigUint {
        &self.exponent_lower_bound
    }

    /// Returns the strategy that generated this accumulation report.
    pub fn strategy(&self) -> GroupExponentStrategy {
        GroupExponentStrategy::RandomPoints {
            max_samples: self.samples_requested,
            point_order_strategy: self.point_order_strategy.clone(),
        }
    }

    /// Returns the recorded accumulation steps.
    pub fn steps(&self) -> &[ExponentAccumulationStep<P>] {
        &self.steps
    }
}
