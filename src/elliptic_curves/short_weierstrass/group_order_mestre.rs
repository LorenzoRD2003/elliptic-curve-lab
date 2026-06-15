use crate::elliptic_curves::frobenius::{
    FrobeniusTrace, GroupOrderReport, MestreConfig, MestreGroupOrderReport, MestreSide,
    MestreStepReport,
};
use crate::elliptic_curves::isomorphisms::TwistKind;
use crate::elliptic_curves::traits::{EnumerableCurveModel, PointIndexSampler};
use crate::elliptic_curves::{
    CurveError, HasseInterval, ShortWeierstrassCurve, ShortWeierstrassQuadraticTwist,
};
use crate::fields::{
    EnumerableFiniteField, FiniteField, FiniteFieldDescriptor, QuadraticCharacterFiniteField,
    SqrtField,
};
use crate::numerics::{NormalizedPrimePowerFactorization, integer_arithmetic::lcm_biguint};

use num_bigint::BigUint;
use std::hash::Hash;

struct MestreFinalizeInput {
    config: MestreConfig,
    base_field: FiniteFieldDescriptor,
    prime: u128,
    interval: crate::elliptic_curves::HasseInterval,
    original_lower_bound: BigUint,
    twist_lower_bound: BigUint,
    steps: Vec<MestreStepReport>,
}

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
        let interval = HasseInterval::for_q(prime)?;
        let twist_curve = self.select_genuine_quadratic_twist_for_mestre()?;

        let mut original_lower_bound = BigUint::from(1u8);
        let mut twist_lower_bound = BigUint::from(1u8);
        let mut next_side = MestreSide::Original;
        let mut steps = Vec::new();
        let mut iterations = 0usize;

        while self
            .unique_multiple_for_side(
                &interval,
                MestreSide::Original,
                &original_lower_bound,
                &twist_lower_bound,
            )
            .is_none()
            && self
                .unique_multiple_for_side(
                    &interval,
                    MestreSide::QuadraticTwist,
                    &original_lower_bound,
                    &twist_lower_bound,
                )
                .is_none()
        {
            if config
                .max_iterations()
                .is_some_and(|max_iterations| iterations >= max_iterations)
            {
                return Err(CurveError::MestreIterationCapReached {
                    max_iterations: config
                        .max_iterations()
                        .expect("checked presence above for the capped branch"),
                });
            }

            let curve = self.curve_for_side(next_side, &twist_curve);
            let lower_bound = self.lower_bound_for_side_mut(
                next_side,
                &mut original_lower_bound,
                &mut twist_lower_bound,
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
            steps.push(MestreStepReport::new(
                next_side,
                annihilating_multiple,
                point_order_report,
                lower_bound.clone(),
            ));

            next_side = self.next_mestre_side(next_side);
            iterations += 1;
        }

        self.finalize_mestre_group_order_report(MestreFinalizeInput {
            config,
            base_field,
            prime,
            interval,
            original_lower_bound,
            twist_lower_bound,
            steps,
        })
    }

    fn unique_mestre_multiple(
        &self,
        interval: &HasseInterval,
        lower_bound: &BigUint,
    ) -> Option<u128> {
        lower_bound
            .try_into()
            .ok()
            .and_then(|bound: u128| interval.unique_multiple_of(bound))
    }

    fn curve_for_side<'a>(
        &'a self,
        side: MestreSide,
        twist_curve: &'a ShortWeierstrassCurve<F>,
    ) -> &'a ShortWeierstrassCurve<F> {
        match side {
            MestreSide::Original => self,
            MestreSide::QuadraticTwist => twist_curve,
        }
    }

    fn lower_bound_for_side_mut<'a>(
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

    fn unique_multiple_for_side(
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

    fn next_mestre_side(&self, side: MestreSide) -> MestreSide {
        match side {
            MestreSide::Original => MestreSide::QuadraticTwist,
            MestreSide::QuadraticTwist => MestreSide::Original,
        }
    }

    fn validate_mestre_prime_field(&self) -> Result<(FiniteFieldDescriptor, u128), CurveError> {
        if F::extension_degree().get() != 1 {
            return Err(CurveError::MestreRequiresPrimeField {
                extension_degree: F::extension_degree().get(),
            });
        }

        if F::characteristic() <= 229 {
            return Err(CurveError::MestrePrimeTooSmall {
                characteristic: F::characteristic(),
            });
        }

        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: F::characteristic(),
                extension_degree: F::extension_degree().get(),
            })?;
        let prime =
            base_field
                .cardinality()
                .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                    characteristic: F::characteristic(),
                    extension_degree: F::extension_degree().get(),
                })?;

        Ok((base_field, prime))
    }

    fn select_genuine_quadratic_twist_for_mestre(
        &self,
    ) -> Result<ShortWeierstrassCurve<F>, CurveError> {
        for candidate in F::elements() {
            if F::is_zero(&candidate) {
                continue;
            }

            let Ok(package) = ShortWeierstrassQuadraticTwist::new(self.clone(), candidate) else {
                continue;
            };
            if package.kind() == TwistKind::Quadratic {
                return Ok(package.twist().clone());
            }
        }

        Err(CurveError::MestreQuadraticTwistUnavailable)
    }

    fn finalize_mestre_group_order_report(
        &self,
        input: MestreFinalizeInput,
    ) -> Result<GroupOrderReport, CurveError> {
        let (resolved_side, original_order, twist_order) = self.resolve_mestre_group_orders(
            &input.interval,
            input.prime,
            &input.original_lower_bound,
            &input.twist_lower_bound,
        );

        let original_trace = FrobeniusTrace::from_order(
            input.base_field.clone(),
            u64::try_from(original_order).expect("Mestre over u64 prime fields should fit in u64"),
        )?;
        let twist_trace = FrobeniusTrace::from_order(
            input.base_field,
            u64::try_from(twist_order).expect("Mestre over u64 prime fields should fit in u64"),
        )?;

        Ok(GroupOrderReport::MestreFp(Box::new(
            MestreGroupOrderReport::new(
                input.config,
                resolved_side,
                original_trace,
                twist_trace,
                input.steps,
            ),
        )))
    }

    fn resolve_mestre_group_orders(
        &self,
        interval: &HasseInterval,
        prime: u128,
        original_lower_bound: &BigUint,
        twist_lower_bound: &BigUint,
    ) -> (MestreSide, u128, u128) {
        if let Some(original_order) = self.unique_multiple_for_side(
            interval,
            MestreSide::Original,
            original_lower_bound,
            twist_lower_bound,
        ) {
            (
                MestreSide::Original,
                original_order,
                self.complementary_group_order_from_twist_sum(prime, original_order),
            )
        } else {
            let twist_order = self
                .unique_multiple_for_side(
                    interval,
                    MestreSide::QuadraticTwist,
                    original_lower_bound,
                    twist_lower_bound,
                )
                .expect("the loop exits only when one Mestre side is unique");
            (
                MestreSide::QuadraticTwist,
                self.complementary_group_order_from_twist_sum(prime, twist_order),
                twist_order,
            )
        }
    }

    fn complementary_group_order_from_twist_sum(&self, prime: u128, known_order: u128) -> u128 {
        prime
            .checked_mul(2)
            .and_then(|double_prime| double_prime.checked_add(2))
            .and_then(|order_sum| order_sum.checked_sub(known_order))
            .expect("2p + 2 - known twist-related order should stay in range")
    }
}
