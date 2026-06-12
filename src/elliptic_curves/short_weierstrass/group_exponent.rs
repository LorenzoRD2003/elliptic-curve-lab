use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    affine::AffinePoint,
    frobenius::{PointCountReport, PointCountStrategy},
    short_weierstrass::PointOrderReport,
    traits::{CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, PointIndexSampler},
};
use crate::fields::{EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField};
use crate::numerics::integer_arithmetic::lcm_biguint;
use num_bigint::BigUint;
use num_traits::ToPrimitive;

use super::PointOrderStrategy;

/// Public strategy choices for recovering or estimating `λ(E(F_q))`.
///
/// For a finite abelian group `G`, the exponent `λ(G) = lcm({|g| : g ∈ G})`
/// is also the maximum element order.
///
/// The current implementation distinguishes:
/// - [`Self::Exhaustive`], which computes the exact exponent on a tiny
///   enumerable curve group
/// - [`Self::RandomPoints`], which samples points with replacement, computes
///   their exact orders through one requested point-order strategy, and
///   accumulates the running least common multiple as a candidate for `λ(G)`
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupExponentStrategy {
    Exhaustive,
    RandomPoints {
        max_samples: usize,
        point_order_strategy: PointOrderStrategy,
    },
}

/// One random-point step in the running `lcm` accumulation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentAccumulationStep<P> {
    point: P,
    point_order_report: PointOrderReport<P>,
    accumulated_lcm: BigUint,
}

impl<P> ExponentAccumulationStep<P> {
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
    fn from_steps(
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

    /// Returns the recorded accumulation steps.
    pub fn steps(&self) -> &[ExponentAccumulationStep<P>] {
        &self.steps
    }
}

/// Shared group-exponent report returned by the unified curve-side API.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupExponentReport<P> {
    Exhaustive(BigUint),
    RandomPoints(Box<ExponentAccumulationReport<P>>),
}

impl<P> GroupExponentReport<P> {
    /// Returns the strategy used to build this report.
    pub fn strategy(&self) -> GroupExponentStrategy {
        match self {
            Self::Exhaustive(_) => GroupExponentStrategy::Exhaustive,
            Self::RandomPoints(report) => GroupExponentStrategy::RandomPoints {
                max_samples: report.samples_requested(),
                point_order_strategy: report.point_order_strategy().clone(),
            },
        }
    }

    /// Returns the best current lower bound for the exponent.
    ///
    /// For the exhaustive route this is exact. For the random-point route this
    /// is the running `lcm` lower bound accumulated from sampled point orders.
    pub fn exponent_lower_bound(&self) -> &BigUint {
        match self {
            Self::Exhaustive(exponent) => exponent,
            Self::RandomPoints(report) => report.exponent_lower_bound(),
        }
    }

    /// Returns the exact exponent when the chosen route computes it
    /// exhaustively.
    pub fn exact_exponent(&self) -> Option<&BigUint> {
        match self {
            Self::Exhaustive(exponent) => Some(exponent),
            Self::RandomPoints(_) => None,
        }
    }
}

/// Point-count-side verification of one accumulated exponent lower bound.
///
/// This report does not certify the exponent itself. It records whether the
/// Hasse interval attached to one chosen point-count route contains a unique
/// multiple of the supplied lower bound, which would force one group order
/// `#E(F_q)` compatible with that lower bound.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentLowerBoundPointCountVerification {
    exponent_lower_bound: BigUint,
    point_count: PointCountReport,
}

impl ExponentLowerBoundPointCountVerification {
    fn new(exponent_lower_bound: BigUint, point_count: PointCountReport) -> Self {
        Self {
            exponent_lower_bound,
            point_count,
        }
    }

    /// Returns the lower bound being checked.
    pub fn exponent_lower_bound(&self) -> &BigUint {
        &self.exponent_lower_bound
    }

    /// Returns the point-count report that supplied the Hasse interval.
    pub fn point_count(&self) -> &PointCountReport {
        &self.point_count
    }

    /// Returns the unique multiple of the lower bound in `H(q)`, if one
    /// exists.
    ///
    /// When this is `Some(N)`, the Hasse interval for the chosen point-count
    /// route contains exactly one multiple of the lower bound, namely `N`.
    /// Since the point-count report already knows the true curve order, this
    /// is best read as a consistency-and-uniqueness witness for `#E(F_q)`,
    /// not as a certification that the exponent itself equals `N`.
    fn unique_group_order_multiple_in_hasse_interval(&self) -> Option<u128> {
        self.exponent_lower_bound.to_u128().and_then(|lower_bound| {
            self.point_count
                .hasse_interval()
                .unique_multiple_of(lower_bound)
        })
    }

    /// Returns the verified group order when the Hasse interval forces one
    /// unique multiple of the lower bound.
    pub fn verified_group_order(&self) -> Option<u128> {
        self.unique_group_order_multiple_in_hasse_interval()
    }
}

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
                    steps.push(ExponentAccumulationStep {
                        point,
                        point_order_report,
                        accumulated_lcm: accumulated_lcm.clone(),
                    });
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

    /// Verifies one accumulated exponent lower bound against a chosen
    /// point-count route.
    ///
    /// This method is intentionally separate from [`Self::group_exponent_by`]:
    /// the random-point exponent route stays a pure lower-bound accumulator,
    /// while this helper uses one explicit [`PointCountStrategy`] to ask
    /// whether the resulting Hasse interval `H(q)` contains a unique multiple
    /// of that lower bound.
    ///
    /// If the returned report has `verified_group_order = Some(N)`, then the
    /// Hasse interval for the chosen point-count route contains exactly one
    /// multiple of the lower bound, namely `N`. This certifies one possible
    /// group order `#E(F_q)`, not the exponent itself.
    ///
    /// The intended workflow is to pass an
    /// [`ExponentAccumulationReport<AffinePoint<F>>`] produced from this same
    /// curve. The method rejects obviously incompatible reports whose sampled
    /// points do not lie on the current curve.
    pub fn verify_exponent_lower_bound_by_point_count(
        &self,
        accumulation: &ExponentAccumulationReport<AffinePoint<F>>,
        strategy: PointCountStrategy,
    ) -> Result<ExponentLowerBoundPointCountVerification, CurveError> {
        for step in accumulation.steps() {
            if !self.contains(step.point()) {
                return Err(CurveError::PointNotOnCurve);
            }
        }

        Ok(ExponentLowerBoundPointCountVerification::new(
            accumulation.exponent_lower_bound().clone(),
            self.count_points(strategy)?,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ExponentLowerBoundPointCountVerification, GroupExponentReport, GroupExponentStrategy,
    };
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, EnumerableCurveModel, PointCountStrategy, PointOrderReport,
        PointOrderStrategy, ShortWeierstrassCurve,
    };
    use crate::fields::{
        EnumerableFiniteField, Field, FiniteField, Fp, QuadraticCharacterFiniteField, SqrtField,
    };
    use num_bigint::BigUint;

    type F7 = Fp<7>;
    type F5 = Fp<5>;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    fn f7_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f5_curve() -> ShortWeierstrassCurve<F5> {
        ShortWeierstrassCurve::<F5>::new(F5::from_i64(0), F5::from_i64(1)).expect("valid curve")
    }

    fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
        f7_curve()
            .point(F7::from_i64(x), F7::from_i64(y))
            .expect("point should lie on the curve")
    }

    fn f5_point(x: i64, y: i64) -> AffinePoint<F5> {
        f5_curve()
            .point(F5::from_i64(x), F5::from_i64(y))
            .expect("point should lie on the curve")
    }

    fn sampler_from_indices(indices: Vec<usize>) -> impl crate::elliptic_curves::PointIndexSampler {
        let mut indices = indices.into_iter();
        move |upper_bound: usize| indices.next().filter(|index| *index < upper_bound)
    }

    fn point_index<F>(curve: &ShortWeierstrassCurve<F>, point: &AffinePoint<F>) -> usize
    where
        F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField,
    {
        curve
            .points()
            .iter()
            .position(|candidate| candidate == point)
            .expect("sample point should appear in the enumerated group")
    }

    #[test]
    fn report_strategy_preserves_the_selected_route() {
        let curve = f7_curve();
        let mut sampler = |_| Some(0usize);

        let exhaustive_report = curve
            .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
            .expect("exhaustive exponent route should succeed");

        assert_eq!(
            exhaustive_report.strategy(),
            GroupExponentStrategy::Exhaustive
        );

        let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(2, 1))]);
        let random_report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy: PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("random-point route should succeed");

        assert_eq!(
            random_report.strategy(),
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            }
        );
    }

    #[test]
    fn group_exponent_by_exhaustive_recovers_the_exact_small_group_exponent() {
        let curve = f7_curve();
        let mut sampler = |_| Some(0usize);

        let report = curve
            .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
            .expect("exhaustive exponent route should succeed");

        assert_eq!(report.strategy(), GroupExponentStrategy::Exhaustive);
        assert_eq!(report.exponent_lower_bound(), &bu(6));
        assert_eq!(report.exact_exponent(), Some(&bu(6)));

        let GroupExponentReport::Exhaustive(exponent) = report else {
            panic!("expected the exhaustive route to preserve its variant");
        };
        assert_eq!(exponent, bu(6));
    }

    #[test]
    fn group_exponent_by_random_points_accumulates_lcms_of_sampled_point_orders() {
        let curve = f7_curve();
        let point_order_two = f7_point(6, 0);
        let point_order_six = f7_point(2, 1);
        let mut sampler = sampler_from_indices(vec![
            point_index(&curve, &point_order_two),
            point_index(&curve, &point_order_six),
        ]);

        let report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 2,
                    point_order_strategy: PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("random-point accumulation should succeed");

        assert_eq!(
            report.strategy(),
            GroupExponentStrategy::RandomPoints {
                max_samples: 2,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            }
        );
        assert_eq!(report.exponent_lower_bound(), &bu(6));
        assert_eq!(report.exact_exponent(), None);

        let GroupExponentReport::RandomPoints(report) = report else {
            panic!("expected the random-point route to preserve its variant");
        };

        assert_eq!(report.samples_requested(), 2);
        assert_eq!(report.samples_taken(), 2);
        assert!(report.completed_requested_samples());
        assert_eq!(
            report.point_order_strategy(),
            &PointOrderStrategy::Exhaustive
        );
        assert_eq!(report.exponent_lower_bound(), &bu(6));
        assert_eq!(report.steps().len(), 2);
        assert_eq!(report.steps()[0].point(), &point_order_two);
        assert_eq!(report.steps()[0].accumulated_lcm(), &bu(2));
        assert_eq!(report.steps()[1].point(), &point_order_six);
        assert_eq!(report.steps()[1].accumulated_lcm(), &bu(6));
    }

    #[test]
    fn group_exponent_by_random_points_preserves_the_point_order_route() {
        let curve = f7_curve();
        let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(2, 1))]);

        let report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy: PointOrderStrategy::HasseIntervalNaive {
                        point_count_strategy: PointCountStrategy::Auto,
                    },
                },
                &mut sampler,
            )
            .expect("Hasse-driven point-order route should compose into exponent accumulation");

        let GroupExponentReport::RandomPoints(report) = report else {
            panic!("expected the random-point route to preserve its variant");
        };

        match report.steps()[0].point_order_report() {
            PointOrderReport::HasseIntervalNaive(step_report) => {
                assert_eq!(
                    step_report.point_count().strategy(),
                    PointCountStrategy::QuadraticCharacter
                );
                assert_eq!(step_report.exact_order(), &bu(6));
            }
            other => panic!("expected HasseIntervalNaive point-order report, got {other:?}"),
        }
    }

    #[test]
    fn group_exponent_by_random_points_reports_early_sampler_exhaustion_honestly() {
        let curve = f7_curve();
        let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(6, 0))]);

        let report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 3,
                    point_order_strategy: PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("short sampling run should still build an honest report");

        let GroupExponentReport::RandomPoints(report) = report else {
            panic!("expected the random-point route to preserve its variant");
        };

        assert_eq!(report.samples_requested(), 3);
        assert_eq!(report.samples_taken(), 1);
        assert!(!report.completed_requested_samples());
        assert_eq!(report.exponent_lower_bound(), &bu(2));
        assert_eq!(report.steps().len(), 1);
    }

    #[test]
    fn point_count_verification_reports_when_the_hasse_interval_is_still_ambiguous() {
        let curve = f7_curve();
        let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(2, 1))]);

        let report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy: PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("accumulation should succeed");

        let GroupExponentReport::RandomPoints(accumulation) = report else {
            panic!("expected the random-point route to preserve its variant");
        };

        let verification = curve
            .verify_exponent_lower_bound_by_point_count(&accumulation, PointCountStrategy::Auto)
            .expect("point-count-side verification should succeed");

        assert_eq!(
            verification,
            ExponentLowerBoundPointCountVerification::new(
                bu(6),
                curve
                    .count_points(PointCountStrategy::Auto)
                    .expect("point count should succeed"),
            )
        );
        assert_eq!(verification.verified_group_order(), None);
    }

    #[test]
    fn point_count_verification_can_force_one_unique_group_order_in_hasse_interval() {
        let curve = f5_curve();
        let mut sampler = sampler_from_indices(vec![point_index(&curve, &f5_point(2, 2))]);

        let report = curve
            .group_exponent_by(
                GroupExponentStrategy::RandomPoints {
                    max_samples: 1,
                    point_order_strategy: PointOrderStrategy::Exhaustive,
                },
                &mut sampler,
            )
            .expect("accumulation should succeed");

        let GroupExponentReport::RandomPoints(accumulation) = report else {
            panic!("expected the random-point route to preserve its variant");
        };

        let verification = curve
            .verify_exponent_lower_bound_by_point_count(
                &accumulation,
                PointCountStrategy::Exhaustive,
            )
            .expect("point-count-side verification should succeed");

        assert_eq!(verification.exponent_lower_bound(), &bu(6));
        assert_eq!(verification.point_count().curve_order(), 6);
        assert_eq!(verification.verified_group_order(), Some(6));
        assert_eq!(
            verification.unique_group_order_multiple_in_hasse_interval(),
            Some(6)
        );
    }
}
