use num_bigint::BigUint;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::point_order::{
        ExhaustivePointOrderReport, HasseIntervalPointOrderReport, PointOrderReport,
        PointOrderStrategy,
    },
    traits::{CurveModel, FiniteGroupCurveModel, HasseIntervalSearchCurveModel},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};
use crate::numerics::NormalizedPrimePowerFactorization;

impl<F: FiniteField + EnumerableFiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    /// Recovers the exact order of one point by one requested strategy.
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
                Ok(PointOrderReport::Exhaustive(
                    ExhaustivePointOrderReport::new(exact_order),
                ))
            }
            PointOrderStrategy::FromKnownMultiple {
                multiple,
                factorization,
            } => self
                .point_order_from_multiple(point, multiple, &factorization)
                .map(PointOrderReport::FromKnownMultiple),
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy,
            } => {
                let group_order_report = self.group_order_by_small_field(group_order_strategy)?;
                let multiple_search = self.find_annihilating_multiple_in_interval_naive(
                    point,
                    group_order_report.hasse_interval(),
                )?;
                let Some(multiple) = multiple_search.first_annihilating_multiple() else {
                    return Err(CurveError::NoAnnihilatingMultipleInHasseInterval {
                        lower: multiple_search.interval().lower(),
                        upper: multiple_search.interval().upper(),
                    });
                };

                let multiple_biguint = multiple.clone();
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
                        group_order_report,
                        multiple_search,
                        order_from_multiple,
                    },
                )))
            }
        }
    }
}
