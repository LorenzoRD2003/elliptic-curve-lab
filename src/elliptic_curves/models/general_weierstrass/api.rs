use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        group_order::{
            FiniteFieldGroupOrderStrategy, GroupOrderReport, SmallFieldGroupOrderStrategy,
        },
    },
    group_algorithms::CyclicPrimaryOrderGroupCurveModel,
    short_weierstrass::{
        group_exponent::{
            ExponentAccumulationReport, ExponentAccumulationStep, GroupExponentReport,
            GroupExponentStrategy,
        },
        point_order::{
            ExhaustivePointOrderReport, HasseIntervalPointOrderReport,
            PointOrderFromMultipleReport, PointOrderReductionStep, PointOrderReport,
            PointOrderStrategy,
        },
    },
    traits::{
        BigScalarGroupCurveModel, CurveModel, CurveModelConversion, EnumerableCurveModel,
        FiniteGroupCurveModel, FrobeniusTraceCurveModel, HasseIntervalSearchCurveModel,
        PointIndexSampler,
    },
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};
use crate::numerics::{
    NormalizedPrimePowerFactorization, PrimePowerTable, integer_arithmetic::lcm_biguint,
};

impl<F: FiniteField> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    /// Computes `#E(F_q)` through one finite-field-capable route.
    ///
    /// Current status:
    /// - `Auto` and `Schoof` are delegated to the short-Weierstrass companion
    ///   when the classical reduction exists
    /// - this wrapper is therefore unavailable in characteristics `2` and `3`
    ///   until the general model gains its own non-enumerative finite-field
    ///   counting route
    pub fn group_order_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        let conversion = self.conversion_to_short_weierstrass()?;
        conversion.target().group_order_by(strategy)
    }

    /// Recovers the Frobenius trace through one finite-field-capable
    /// group-order route.
    ///
    /// Current status:
    /// - delegated to the short-Weierstrass companion
    /// - unavailable in characteristics `2` and `3` for now
    pub fn frobenius_trace_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }
}

impl<F> GeneralWeierstrassCurve<F>
where
    F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField,
    F::Elem: Clone,
{
    fn validate_point_order_from_multiple_inputs(
        &self,
        point: &AffinePoint<F>,
        multiple: &BigUint,
    ) -> Result<(), CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        if multiple.is_zero() {
            return Err(CurveError::InvalidPointOrderMultiple {
                multiple: multiple.clone(),
            });
        }

        let image = self.mul_scalar_biguint(point, multiple)?;
        if self.is_identity(&image) {
            Ok(())
        } else {
            Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
                multiple: multiple.clone(),
            })
        }
    }

    fn recover_point_order_from_normalized_factorization(
        &self,
        point: &AffinePoint<F>,
        supplied_multiple: BigUint,
        normalized_factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        let mut remaining_multiple = supplied_multiple.clone();
        let mut exact_order = BigUint::from(1u8);
        let mut steps = Vec::with_capacity(normalized_factorization.len());

        for (prime, exponent_in_multiple) in normalized_factorization {
            let powers = PrimePowerTable::up_through(prime, *exponent_in_multiple);
            let prime_power = powers.power(*exponent_in_multiple);
            let cofactor = &remaining_multiple / prime_power;
            let primary_component = if cofactor == BigUint::from(1u8) {
                point.clone()
            } else {
                self.mul_scalar_biguint(point, &cofactor)?
            };

            let local_report = self.recover_cyclic_primary_order(&primary_component, &powers)?;
            let removed_exponent = local_report.removed_exponent();
            let local_exact_power = powers.power(local_report.exact_exponent());
            let removed_power = powers.power(removed_exponent);

            exact_order *= local_exact_power;
            remaining_multiple /= removed_power;

            steps.push(PointOrderReductionStep::new(
                local_report.prime().clone(),
                local_report.exponent_bound(),
                removed_exponent,
                remaining_multiple.clone(),
            ));
        }

        Ok(PointOrderFromMultipleReport::new(
            supplied_multiple,
            exact_order,
            steps,
        ))
    }

    fn point_order_from_multiple(
        &self,
        point: &AffinePoint<F>,
        multiple: BigUint,
        factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        self.validate_point_order_from_multiple_inputs(point, &multiple)?;
        let normalized_factorization =
            NormalizedPrimePowerFactorization::checked(&multiple, factorization)
                .map_err(|_| CurveError::InvalidPointOrderMultipleFactorization {
                    multiple: multiple.clone(),
                })?
                .into_factors();
        self.recover_point_order_from_normalized_factorization(
            point,
            multiple,
            &normalized_factorization,
        )
    }

    /// Computes `#E(F_q)` through one route that is specific to small
    /// enumerable finite fields.
    ///
    /// Current status:
    /// - `Exhaustive` is native to `GeneralWeierstrassCurve<F>` through direct
    ///   point enumeration
    /// - `Auto`, `QuadraticCharacter`, and `Schoof` reuse the short companion
    ///   when available
    /// - `Auto` falls back to the native exhaustive route in characteristics
    ///   `2` and `3`, where the short reduction is unavailable
    pub fn group_order_by_small_field(
        &self,
        strategy: SmallFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError>
    where
        Self: FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
        match strategy {
            SmallFieldGroupOrderStrategy::Exhaustive => {
                FrobeniusTraceCurveModel::frobenius_trace(self)
                    .map(GroupOrderReport::ExhaustiveTrace)
            }
            SmallFieldGroupOrderStrategy::Auto if matches!(F::characteristic(), 2 | 3) => {
                FrobeniusTraceCurveModel::frobenius_trace(self)
                    .map(GroupOrderReport::ExhaustiveTrace)
            }
            SmallFieldGroupOrderStrategy::Auto
            | SmallFieldGroupOrderStrategy::QuadraticCharacter
            | SmallFieldGroupOrderStrategy::Schoof => {
                let conversion = self.conversion_to_short_weierstrass()?;
                conversion.target().group_order_by_small_field(strategy)
            }
        }
    }

    /// Recovers the exact order of one point by one requested strategy.
    ///
    /// Current status:
    /// - `Exhaustive` is native
    /// - `FromKnownMultiple` is native through the shared big-scalar and
    ///   cyclic-primary helpers
    /// - `HasseIntervalNaive` is native once the chosen small-field
    ///   group-order route produces a Hasse interval
    pub fn point_order_by(
        &self,
        point: &AffinePoint<F>,
        strategy: PointOrderStrategy,
    ) -> Result<PointOrderReport<AffinePoint<F>>, CurveError>
    where
        Self: FiniteGroupCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
            + FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
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

                let multiple_biguint = BigUint::from(multiple);
                let factorization = NormalizedPrimePowerFactorization::factor(&multiple_biguint)
                    .expect("an annihilating multiple in H(q) should admit a prime factorization")
                    .into_factors();
                let order_from_multiple = self.recover_point_order_from_normalized_factorization(
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

    /// Recovers or estimates `λ(E(F_q))` by one requested strategy.
    ///
    /// Current status:
    /// - `Exhaustive` is native
    /// - `RandomPoints` is native once the chosen point-order route is
    ///   available
    pub fn group_exponent_by<S: PointIndexSampler>(
        &self,
        strategy: GroupExponentStrategy,
        sampler: &mut S,
    ) -> Result<GroupExponentReport<AffinePoint<F>>, CurveError>
    where
        Self: FiniteGroupCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
            + FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
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
