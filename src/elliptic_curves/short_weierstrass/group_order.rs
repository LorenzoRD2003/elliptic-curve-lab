use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{CharacterSumPointCount, FrobeniusTrace, GroupOrderReport, GroupOrderStrategy},
    traits::{FrobeniusTraceCurveModel, PointIndexSampler},
};
use crate::fields::{
    EnumerableFiniteField, FiniteField, FiniteFieldDescriptor, QuadraticCharacterFiniteField,
    SqrtField,
};
use std::hash::Hash;

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    fn group_order_by_without_sampler(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        match strategy {
            GroupOrderStrategy::Auto | GroupOrderStrategy::QuadraticCharacter => self
                .group_order_by_quadratic_character()
                .map(GroupOrderReport::QuadraticCharacter),
            GroupOrderStrategy::Exhaustive => FrobeniusTraceCurveModel::frobenius_trace(self)
                .map(GroupOrderReport::ExhaustiveTrace),
            GroupOrderStrategy::MestreFp(_) => Err(CurveError::GroupOrderStrategyRequiresSampler {
                strategy: "MestreFp",
            }),
        }
    }

    /// Computes `#E(F_q)` using one requested public strategy.
    ///
    /// This is the user-facing finite-field group-order entry point for the
    /// current short-Weierstrass model for deterministic routes.
    ///
    /// Complexity:
    /// - [`GroupOrderStrategy::Exhaustive`]: dominated by full rational-point
    ///   enumeration
    /// - [`GroupOrderStrategy::QuadraticCharacter`]: `Θ(q)` right-hand-side
    ///   evaluations and quadratic-character queries
    /// - [`GroupOrderStrategy::Auto`]: currently the same as the
    ///   quadratic-character route
    ///
    /// Strategies that need a point sampler, such as [`GroupOrderStrategy::MestreFp`],
    /// must use [`Self::group_order_by_with_sampler`].
    pub fn group_order_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        self.group_order_by_without_sampler(strategy)
    }

    /// Computes `#E(F_q)` using one requested public strategy, including
    /// sampler-driven routes such as Mestre's prime-field algorithm.
    ///
    /// Complexity:
    /// - deterministic strategies delegate to [`Self::group_order_by`]
    /// - [`GroupOrderStrategy::MestreFp`]: expected `Θ(p^(1/4) M(log p))` bit
    ///   complexity in the current BSGS-backed implementation over `F_p`.
    pub fn group_order_by_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: GroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<GroupOrderReport, CurveError>
    where
        F::Elem: Hash,
    {
        match strategy {
            GroupOrderStrategy::MestreFp(config) => self.group_order_by_mestre_fp(config, sampler),
            other => self.group_order_by_without_sampler(other),
        }
    }

    /// Internal short-Weierstrass-specific `Θ(q)` character-sum count.
    ///
    /// The public entry point for callers is [`Self::group_order_by`]. This
    /// helper stays crate-private so the public API has one primary counting
    /// door while still letting internal tests and visualizations exercise the
    /// specific route directly.
    ///
    /// Formula:
    /// `#E(F_q) = q + 1 + Σ_{x ∈ F_q} χ(x^3 + Ax + B)`.
    ///
    /// Complexity:
    /// `Θ(q)` evaluations of `x^3 + Ax + B` and `Θ(q)` quadratic-character
    /// queries over represented field elements.
    pub(crate) fn group_order_by_quadratic_character(
        &self,
    ) -> Result<CharacterSumPointCount, CurveError> {
        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: F::characteristic(),
                extension_degree: F::extension_degree().get(),
            })?;

        let mut character_sum = 0i128;
        for x in F::elements() {
            let rhs = self.rhs_value(&x);
            let value = F::quadratic_character_of(&rhs).map_err(|_| {
                CurveError::UnsupportedCharacterSumPointCount {
                    characteristic: F::characteristic(),
                    extension_degree: F::extension_degree().get(),
                }
            })?;
            character_sum += value.as_i128();
        }

        CharacterSumPointCount::new(base_field, character_sum)
    }

    /// Recovers the Frobenius trace through one requested group-order
    /// strategy for deterministic routes.
    pub fn frobenius_trace_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }

    /// Recovers the Frobenius trace through one requested group-order
    /// strategy, including sampler-driven routes.
    pub fn frobenius_trace_by_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: GroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<FrobeniusTrace, CurveError>
    where
        F::Elem: Hash,
    {
        self.group_order_by_with_sampler(strategy, sampler)?
            .to_frobenius_trace()
    }
}

#[cfg(test)]
mod tests {
    use super::ShortWeierstrassCurve;
    use crate::elliptic_curves::{
        CurveError, EnumerableCurveModel, FiniteGroupCurveModel, GroupOrderReport,
        GroupOrderStrategy, MestreConfig, MestreSide, ShortWeierstrassQuadraticTwist, TwistKind,
    };
    use crate::fields::{EnumerableFiniteField, Field, Fp};
    use crate::proptest_support::config::CurveStrategyConfig;
    use crate::proptest_support::elliptic_curves::arb_nonsingular_curve;
    use proptest::prelude::*;
    use std::collections::HashMap;

    type F43 = Fp<43>;
    type F241 = Fp<241>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F43Sqrt2MestreSpec,
        field: F43Sqrt2Mestre,
        base: F43,
        non_residue: 2,
        name: "F43(sqrt(2)) for Mestre tests",
    );

    fn f241_curve() -> ShortWeierstrassCurve<F241> {
        ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
            .expect("valid F241 curve")
    }

    fn max_order_point_index<F>(curve: &ShortWeierstrassCurve<F>) -> usize
    where
        F: crate::fields::FiniteField
            + crate::fields::EnumerableFiniteField
            + crate::fields::QuadraticCharacterFiniteField
            + crate::fields::SqrtField,
    {
        curve
            .point_orders()
            .into_iter()
            .enumerate()
            .max_by_key(|(_, (_, order))| *order)
            .map(|(index, _)| index)
            .expect("small finite curve should have at least one point")
    }

    fn genuine_twist_curve(curve: &ShortWeierstrassCurve<F241>) -> ShortWeierstrassCurve<F241> {
        F241::elements()
            .into_iter()
            .find_map(|candidate| {
                if F241::is_zero(&candidate) {
                    return None;
                }
                let Ok(package) = ShortWeierstrassQuadraticTwist::new(curve.clone(), candidate)
                else {
                    return None;
                };
                (package.kind() == TwistKind::Quadratic).then(|| package.twist().clone())
            })
            .expect("a prime-field curve should admit a genuine quadratic twist")
    }

    fn sampler_covering_each_curve_by_index() -> impl FnMut(usize) -> Option<usize> {
        let mut next_index_by_upper_bound = HashMap::<usize, usize>::new();
        move |upper_bound: usize| {
            let next_index = next_index_by_upper_bound.entry(upper_bound).or_insert(0);
            let sampled = *next_index % upper_bound;
            *next_index += 1;
            Some(sampled)
        }
    }

    #[test]
    fn mestre_group_order_by_with_sampler_matches_exhaustive_on_a_prime_field_curve() {
        let curve = f241_curve();
        let twist_curve = genuine_twist_curve(&curve);
        let original_index = max_order_point_index(&curve);
        let twist_index = max_order_point_index(&twist_curve);
        let mut requested = vec![original_index, twist_index].into_iter();
        let mut sampler = move |_upper_bound: usize| requested.next().or(Some(original_index));

        let report = curve
            .group_order_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(8)),
                &mut sampler,
            )
            .expect("Mestre route should recover the group order over F241");
        let exhaustive = curve
            .group_order_by(GroupOrderStrategy::Exhaustive)
            .expect("exhaustive group order should succeed");

        assert_eq!(report.curve_order(), exhaustive.curve_order());
        assert_eq!(report.trace(), exhaustive.trace());

        let GroupOrderReport::MestreFp(mestre_report) = report else {
            panic!("Mestre strategy should preserve its own report variant");
        };

        assert!(matches!(
            mestre_report.resolved_side(),
            MestreSide::Original | MestreSide::QuadraticTwist
        ));
        assert!(!mestre_report.steps().is_empty());
    }

    #[test]
    fn mestre_frobenius_trace_with_sampler_matches_exhaustive_trace() {
        let curve = f241_curve();
        let twist_curve = genuine_twist_curve(&curve);
        let original_index = max_order_point_index(&curve);
        let twist_index = max_order_point_index(&twist_curve);
        let mut requested = vec![original_index, twist_index].into_iter();
        let mut sampler = move |_upper_bound: usize| requested.next().or(Some(original_index));

        let mestre_trace = curve
            .frobenius_trace_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(8)),
                &mut sampler,
            )
            .expect("Mestre route should recover the Frobenius trace");

        assert_eq!(
            mestre_trace,
            curve
                .frobenius_trace_by(GroupOrderStrategy::Exhaustive)
                .expect("exhaustive trace should compute")
        );
    }

    #[test]
    fn deterministic_group_order_api_reports_that_mestre_needs_a_sampler() {
        let curve = f241_curve();

        assert_eq!(
            curve.group_order_by(GroupOrderStrategy::MestreFp(MestreConfig::unbounded())),
            Err(CurveError::GroupOrderStrategyRequiresSampler {
                strategy: "MestreFp"
            })
        );
    }

    #[test]
    fn mestre_route_rejects_prime_fields_below_the_theorem_threshold() {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
        let mut sampler = |_upper_bound: usize| Some(0usize);

        assert_eq!(
            curve.group_order_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
                &mut sampler,
            ),
            Err(CurveError::MestrePrimeTooSmall { characteristic: 43 })
        );
    }

    #[test]
    fn mestre_route_rejects_extension_fields() {
        let curve = ShortWeierstrassCurve::<F43Sqrt2Mestre>::new(
            F43Sqrt2Mestre::one(),
            F43Sqrt2Mestre::one(),
        )
        .expect("valid extension-field curve");
        let mut sampler = |_upper_bound: usize| Some(0usize);

        assert_eq!(
            curve.group_order_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
                &mut sampler,
            ),
            Err(CurveError::MestreRequiresPrimeField {
                extension_degree: 2
            })
        );
    }

    #[test]
    fn mestre_route_reports_iteration_cap_reached_before_sampling() {
        let curve = f241_curve();
        let mut sampler = |_upper_bound: usize| Some(0usize);

        assert_eq!(
            curve.group_order_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(0)),
                &mut sampler,
            ),
            Err(CurveError::MestreIterationCapReached { max_iterations: 0 })
        );
    }

    #[test]
    fn mestre_route_reports_sampler_exhaustion() {
        let curve = f241_curve();
        let mut sampler = |_upper_bound: usize| None::<usize>;

        assert_eq!(
            curve.group_order_by_with_sampler(
                GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
                &mut sampler,
            ),
            Err(CurveError::MestreSamplerExhausted)
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        #[test]
        fn property_mestre_matches_exhaustive_group_order_over_f241(
            curve in arb_nonsingular_curve::<241>(CurveStrategyConfig::default()),
        ) {
            let twist_curve = genuine_twist_curve(&curve);
            let max_iterations = 2 * curve.order().max(twist_curve.order());
            let mut sampler = sampler_covering_each_curve_by_index();

            let mestre = curve
                .group_order_by_with_sampler(
                    GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(max_iterations)),
                    &mut sampler,
                )
                .expect("Mestre should recover the group order over F241 after covering both curves");
            let exhaustive = curve
                .group_order_by(GroupOrderStrategy::Exhaustive)
                .expect("exhaustive group order should succeed over F241");

            prop_assert_eq!(mestre.curve_order(), exhaustive.curve_order());
            prop_assert_eq!(mestre.trace(), exhaustive.trace());

            let GroupOrderReport::MestreFp(report) = mestre else {
                panic!("Mestre strategy should preserve its own report variant");
            };

            prop_assert!(matches!(
                report.resolved_side(),
                MestreSide::Original | MestreSide::QuadraticTwist
            ));
            prop_assert!(!report.steps().is_empty());
            prop_assert!(report.steps().len() <= max_iterations);
        }
    }
}
