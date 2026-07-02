use num_bigint::BigUint;
use num_traits::Zero;

use crate::elliptic_curves::{
    CurveError,
    frobenius::group_order::{GroupOrderReport, SmallFieldGroupOrderStrategy},
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
        BigScalarGroupCurveModel, EnumerableCurveModel, FiniteGroupCurveModel,
        HasseIntervalSearchCurveModel, PointIndexSampler,
    },
};
use crate::fields::traits::EnumerableFiniteField;
use crate::numerics::{
    NormalizedPrimePowerFactorization, PrimePowerTable, integer_arithmetic::lcm_biguint,
};

pub(crate) fn validate_point_order_from_multiple_inputs<C: BigScalarGroupCurveModel>(
    curve: &C,
    point: &C::Point,
    multiple: &BigUint,
) -> Result<(), CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }
    if multiple.is_zero() {
        return Err(CurveError::InvalidPointOrderMultiple {
            multiple: multiple.clone(),
        });
    }

    let image = curve.mul_scalar_biguint(point, multiple)?;
    if curve.is_identity(&image) {
        Ok(())
    } else {
        Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
            multiple: multiple.clone(),
        })
    }
}

pub(crate) fn recover_point_order_from_normalized_factorization<
    C: CyclicPrimaryOrderGroupCurveModel,
>(
    curve: &C,
    point: &C::Point,
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
            curve.mul_scalar_biguint(point, &cofactor)?
        };

        let local_report = curve.recover_cyclic_primary_order(&primary_component, &powers)?;
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

pub(crate) fn point_order_from_multiple<C: CyclicPrimaryOrderGroupCurveModel>(
    curve: &C,
    point: &C::Point,
    multiple: BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<PointOrderFromMultipleReport, CurveError> {
    validate_point_order_from_multiple_inputs(curve, point, &multiple)?;
    let normalized_factorization =
        NormalizedPrimePowerFactorization::checked(&multiple, factorization)
            .map_err(|_| CurveError::InvalidPointOrderMultipleFactorization {
                multiple: multiple.clone(),
            })?
            .into_factors();
    recover_point_order_from_normalized_factorization(
        curve,
        point,
        multiple,
        &normalized_factorization,
    )
}

pub(crate) fn point_order_by<C, G>(
    curve: &C,
    point: &C::Point,
    strategy: PointOrderStrategy,
    group_order_by_small_field: G,
) -> Result<PointOrderReport<C::Point>, CurveError>
where
    C: FiniteGroupCurveModel + HasseIntervalSearchCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
    G: Fn(SmallFieldGroupOrderStrategy) -> Result<GroupOrderReport, CurveError>,
{
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    match strategy {
        PointOrderStrategy::Exhaustive => {
            let exact_order = curve
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
        } => point_order_from_multiple(curve, point, multiple, &factorization)
            .map(PointOrderReport::FromKnownMultiple),
        PointOrderStrategy::HasseIntervalNaive {
            group_order_strategy,
        } => {
            let group_order_report = group_order_by_small_field(group_order_strategy)?;
            let multiple_search = curve.find_annihilating_multiple_in_interval_naive(
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
            let order_from_multiple = recover_point_order_from_normalized_factorization(
                curve,
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

pub(crate) fn group_exponent_by<C, S, P>(
    curve: &C,
    strategy: GroupExponentStrategy,
    sampler: &mut S,
    point_order_by: P,
) -> Result<GroupExponentReport<C::Point>, CurveError>
where
    C: FiniteGroupCurveModel + EnumerableCurveModel,
    C::BaseField: EnumerableFiniteField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
    S: PointIndexSampler,
    P: Fn(&C::Point, PointOrderStrategy) -> Result<PointOrderReport<C::Point>, CurveError>,
{
    match strategy {
        GroupExponentStrategy::Exhaustive => Ok(GroupExponentReport::Exhaustive(BigUint::from(
            curve.exponent(),
        ))),
        GroupExponentStrategy::RandomPoints {
            max_samples,
            point_order_strategy,
        } => {
            let mut steps = Vec::with_capacity(max_samples);
            let mut accumulated_lcm = BigUint::from(1u8);

            for _ in 0..max_samples {
                let Some(point) = curve.random_point(sampler) else {
                    break;
                };
                let point_order_report = point_order_by(&point, point_order_strategy.clone())?;
                accumulated_lcm = lcm_biguint(&accumulated_lcm, point_order_report.exact_order());
                steps.push(ExponentAccumulationStep::new(
                    point,
                    point_order_report,
                    accumulated_lcm.clone(),
                ));
            }

            Ok(GroupExponentReport::RandomPoints(Box::new(
                ExponentAccumulationReport::from_steps(max_samples, point_order_strategy, steps),
            )))
        }
    }
}
