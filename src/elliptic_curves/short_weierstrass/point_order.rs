use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    affine::AffinePoint,
    frobenius::{HasseMultipleSearchReport, PointCountReport, PointCountStrategy},
    short_weierstrass::PointOrderFromMultipleReport,
    traits::{CurveModel, FiniteGroupCurveModel, HasseMultipleSearchCurveModel},
};
use crate::fields::{EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField};
use crate::numerics::NormalizedPrimePowerFactorization;
use num_bigint::BigUint;

/// Public strategy choices for recovering the exact order of one point.
///
/// The current educational implementation distinguishes:
///
/// - [`Self::Exhaustive`], which traverses the small ambient finite group
/// - [`Self::FromKnownMultiple`], which starts from one supplied annihilating
///   multiple and peels prime powers
/// - [`Self::HasseIntervalNaive`], which first counts points by one requested
///   public point-count strategy, derives `H(q)` from that report, finds one
///   annihilating multiple in `H(q)`, and then reuses the prime-peeling route
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointOrderStrategy {
    Exhaustive,
    FromKnownMultiple {
        multiple: BigUint,
        factorization: Vec<(BigUint, u32)>,
    },
    HasseIntervalNaive {
        point_count_strategy: PointCountStrategy,
    },
}

/// Route labels for point-order reports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointOrderStrategyKind {
    Exhaustive,
    FromKnownMultiple,
    HasseIntervalNaive,
}

impl PointOrderStrategy {
    /// Returns the route label without the strategy payload.
    pub fn kind(&self) -> PointOrderStrategyKind {
        match self {
            Self::Exhaustive => PointOrderStrategyKind::Exhaustive,
            Self::FromKnownMultiple { .. } => PointOrderStrategyKind::FromKnownMultiple,
            Self::HasseIntervalNaive { .. } => PointOrderStrategyKind::HasseIntervalNaive,
        }
    }
}

/// Report for the exhaustive small-group point-order route.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExhaustivePointOrderReport {
    exact_order: BigUint,
}

impl ExhaustivePointOrderReport {
    /// Returns the recovered exact order.
    pub fn exact_order(&self) -> &BigUint {
        &self.exact_order
    }
}

/// Report for the Hasse-interval route to point order.
///
/// This route first finds one annihilating multiple `M ∈ H(q)` and then
/// recovers `ord(P)` from that `M`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseIntervalPointOrderReport<P> {
    point_count: PointCountReport,
    multiple_search: HasseMultipleSearchReport<P>,
    order_from_multiple: PointOrderFromMultipleReport,
}

impl<P> HasseIntervalPointOrderReport<P> {
    /// Returns the point-count report used to derive `H(q)`.
    pub fn point_count(&self) -> &PointCountReport {
        &self.point_count
    }

    /// Returns the Hasse-interval multiple search report.
    pub fn multiple_search(&self) -> &HasseMultipleSearchReport<P> {
        &self.multiple_search
    }

    /// Returns the follow-up prime-peeling order report.
    pub fn order_from_multiple(&self) -> &PointOrderFromMultipleReport {
        &self.order_from_multiple
    }

    /// Returns the recovered exact order.
    pub fn exact_order(&self) -> &BigUint {
        self.order_from_multiple.exact_order()
    }
}

/// Shared point-order report returned by the unified curve-side order API.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointOrderReport<P> {
    Exhaustive(ExhaustivePointOrderReport),
    FromKnownMultiple(PointOrderFromMultipleReport),
    HasseIntervalNaive(Box<HasseIntervalPointOrderReport<P>>),
}

impl<P> PointOrderReport<P> {
    /// Returns the strategy kind used to build this report.
    pub fn strategy_kind(&self) -> PointOrderStrategyKind {
        match self {
            Self::Exhaustive(_) => PointOrderStrategyKind::Exhaustive,
            Self::FromKnownMultiple(_) => PointOrderStrategyKind::FromKnownMultiple,
            Self::HasseIntervalNaive(_) => PointOrderStrategyKind::HasseIntervalNaive,
        }
    }

    /// Returns the recovered exact order.
    pub fn exact_order(&self) -> &BigUint {
        match self {
            Self::Exhaustive(report) => report.exact_order(),
            Self::FromKnownMultiple(report) => report.exact_order(),
            Self::HasseIntervalNaive(report) => report.exact_order(),
        }
    }
}

impl<F> ShortWeierstrassCurve<F>
where
    F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField,
{
    /// Recovers the exact order of one point by one requested strategy.
    ///
    /// Complexity:
    /// - [`PointOrderStrategy::Exhaustive`]: `Θ(#E(F_q))` group additions in
    ///   the current direct traversal implementation
    /// - [`PointOrderStrategy::FromKnownMultiple`]: delegated to
    ///   [`Self::point_order_from_multiple`]
    /// - [`PointOrderStrategy::HasseIntervalNaive`]: one naive Hasse-interval
    ///   annihilating-multiple search, then one prime-peeling recovery from
    ///   the found multiple
    pub fn point_order_by(
        &self,
        point: &AffinePoint<F>,
        strategy: PointOrderStrategy,
    ) -> Result<PointOrderReport<AffinePoint<F>>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match strategy {
            PointOrderStrategy::Exhaustive => {
                let exact_order = self
                    .point_order(point)
                    .map(BigUint::from)
                    .expect("validated small finite curve points should have an exact order");
                Ok(PointOrderReport::Exhaustive(ExhaustivePointOrderReport {
                    exact_order,
                }))
            }
            PointOrderStrategy::FromKnownMultiple {
                multiple,
                factorization,
            } => self
                .point_order_from_multiple(point, multiple, &factorization)
                .map(PointOrderReport::FromKnownMultiple),
            PointOrderStrategy::HasseIntervalNaive {
                point_count_strategy,
            } => {
                let point_count = self.count_points(point_count_strategy)?;
                let multiple_search = self.find_annihilating_multiple_in_interval_naive(
                    point,
                    point_count.hasse_interval(),
                )?;
                let Some(multiple) = multiple_search.first_annihilating_multiple() else {
                    return Err(CurveError::NoAnnihilatingMultipleInHasseInterval {
                        lower: multiple_search.interval().lower(),
                        upper: multiple_search.interval().upper(),
                    });
                };

                let multiple_biguint = BigUint::from(multiple);
                let factorization = NormalizedPrimePowerFactorization::factor(&multiple_biguint)
                    .expect("an annihilating multiple in H(q) should admit a prime factorization")
                    .into_factors();
                let order_from_multiple = self
                    .point_order_from_multiple_with_trusted_factorization(
                        point,
                        multiple_biguint,
                        &factorization,
                    )?;

                Ok(PointOrderReport::HasseIntervalNaive(Box::new(
                    HasseIntervalPointOrderReport {
                        point_count,
                        multiple_search,
                        order_from_multiple,
                    },
                )))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ExhaustivePointOrderReport, HasseIntervalPointOrderReport, PointOrderReport,
        PointOrderStrategy, PointOrderStrategyKind,
    };
    use crate::elliptic_curves::{AffineCurveModel, PointCountStrategy, ShortWeierstrassCurve};
    use crate::fields::{Field, Fp};
    use crate::numerics::NormalizedPrimePowerFactorization;
    use num_bigint::BigUint;

    type F7 = Fp<7>;

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
    }

    fn f7_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn strategy_kind_preserves_the_selected_route() {
        assert_eq!(
            PointOrderStrategy::Exhaustive.kind(),
            PointOrderStrategyKind::Exhaustive
        );
        assert_eq!(
            PointOrderStrategy::FromKnownMultiple {
                multiple: bu(6),
                factorization: vec![(bu(2), 1), (bu(3), 1)],
            }
            .kind(),
            PointOrderStrategyKind::FromKnownMultiple
        );
        assert_eq!(
            PointOrderStrategy::HasseIntervalNaive {
                point_count_strategy: PointCountStrategy::Auto,
            }
            .kind(),
            PointOrderStrategyKind::HasseIntervalNaive
        );
    }

    #[test]
    fn point_order_by_exhaustive_recovers_the_small_exact_order() {
        let curve = f7_curve();
        let point = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("sample point should lie on the curve");

        let report = curve
            .point_order_by(&point, PointOrderStrategy::Exhaustive)
            .expect("exhaustive route should succeed");

        assert_eq!(report.strategy_kind(), PointOrderStrategyKind::Exhaustive);
        assert_eq!(report.exact_order(), &bu(6));
        assert_eq!(
            report,
            PointOrderReport::Exhaustive(ExhaustivePointOrderReport { exact_order: bu(6) })
        );
    }

    #[test]
    fn point_order_by_from_known_multiple_reuses_the_prime_peeling_report() {
        let curve = f7_curve();
        let point = curve
            .point(F7::from_i64(6), F7::from_i64(0))
            .expect("sample point should lie on the curve");

        let report = curve
            .point_order_by(
                &point,
                PointOrderStrategy::FromKnownMultiple {
                    multiple: bu(6),
                    factorization: vec![(bu(2), 1), (bu(3), 1)],
                },
            )
            .expect("known-multiple route should succeed");

        assert_eq!(
            report.strategy_kind(),
            PointOrderStrategyKind::FromKnownMultiple
        );
        assert_eq!(report.exact_order(), &bu(2));
    }

    #[test]
    fn point_order_by_hasse_interval_naive_composes_search_and_prime_peeling() {
        let curve = f7_curve();
        let point = curve
            .point(F7::from_i64(2), F7::from_i64(1))
            .expect("sample point should lie on the curve");

        let report = curve
            .point_order_by(
                &point,
                PointOrderStrategy::HasseIntervalNaive {
                    point_count_strategy: PointCountStrategy::Auto,
                },
            )
            .expect("Hasse-interval route should succeed");

        assert_eq!(
            report.strategy_kind(),
            PointOrderStrategyKind::HasseIntervalNaive
        );
        assert_eq!(report.exact_order(), &bu(6));

        let PointOrderReport::HasseIntervalNaive(report) = report else {
            panic!("expected the Hasse-interval route to preserve its report variant");
        };

        let HasseIntervalPointOrderReport {
            point_count,
            multiple_search,
            order_from_multiple,
        } = *report;

        assert_eq!(multiple_search.q(), 7);
        assert_eq!(multiple_search.interval(), &point_count.hasse_interval());
        assert_eq!(
            point_count.strategy(),
            PointCountStrategy::QuadraticCharacter
        );
        assert_eq!(multiple_search.first_annihilating_multiple(), Some(6));
        assert_eq!(order_from_multiple.supplied_multiple(), &bu(6));
        assert_eq!(order_from_multiple.exact_order(), &bu(6));
    }

    #[test]
    fn point_order_by_rejects_off_curve_inputs() {
        let curve = f7_curve();
        let invalid =
            crate::elliptic_curves::AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert_eq!(
            curve.point_order_by(&invalid, PointOrderStrategy::Exhaustive),
            Err(crate::elliptic_curves::CurveError::PointNotOnCurve)
        );
    }

    #[test]
    fn hasse_route_uses_the_shared_normalized_prime_power_factorization_surface() {
        assert_eq!(
            NormalizedPrimePowerFactorization::factor(&bu(72))
                .expect("72 should factor")
                .into_factors(),
            vec![(bu(2), 3), (bu(3), 2)]
        );
    }
}
