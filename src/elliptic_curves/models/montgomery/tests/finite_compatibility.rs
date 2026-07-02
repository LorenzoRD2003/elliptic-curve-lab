use crate::fields::traits::*;
use num_bigint::{BigInt, BigUint};
use proptest::prelude::*;

use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        group_order::{
            FiniteFieldGroupOrderStrategy, GroupOrderReport, GroupOrderRoute,
            SmallFieldGroupOrderStrategy,
        },
    },
    short_weierstrass::point_order::{PointOrderReport, PointOrderStrategy},
    traits::{
        AffineCurveModel, CurveModelConversion, EnumerableCurveModel, FiniteAbelianGroupStructure,
        FiniteGroupCurveModel, FrobeniusTraceCurveModel, GroupCurveModel,
    },
};
use crate::proptest_support::{
    config::CurveStrategyConfig,
    elliptic_curves::{arb_montgomery_curve_and_point, arb_nonsingular_montgomery_curve},
};

type F3 = crate::fields::Fp3;
type F5 = crate::fields::Fp5;

use super::shared::{f3_curve, f5_curve};

#[test]
fn finite_group_curve_model_reports_the_expected_structure_in_characteristic_three() {
    let curve = f3_curve();
    let structure = curve.group_structure();

    assert_eq!(curve.order(), 4);
    assert_eq!(curve.exponent(), 4);
    assert_eq!(curve.order_distribution().get(&1), Some(&1));
    assert_eq!(curve.order_distribution().get(&2), Some(&1));
    assert_eq!(curve.order_distribution().get(&4), Some(&2));
    assert_eq!(
        structure,
        FiniteAbelianGroupStructure {
            order: 4,
            exponent: 4,
            cyclic: true,
            invariant_factors: None,
        }
    );
}

#[test]
fn frobenius_trace_curve_model_recovers_the_expected_trace_in_characteristic_three() {
    let curve = f3_curve();
    let trace = curve
        .frobenius_trace()
        .expect("small finite curve should yield a Frobenius trace");

    assert_eq!(trace.curve_order(), BigUint::from(4u8));
    assert_eq!(trace.field_order(), BigUint::from(3u8));
    assert_eq!(trace.trace(), BigInt::from(0));
    assert!(
        curve
            .verify_hasse_bound()
            .expect("Hasse bound verification should succeed")
            .holds()
    );
}

#[test]
fn hasse_search_finds_an_annihilating_multiple_for_a_montgomery_curve_point() {
    let curve = f5_curve();
    let point = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the curve");

    let report = curve
        .find_annihilating_multiple_in_hasse_interval_naive(&point)
        .expect("naive Hasse search should succeed");
    let annihilating_multiple = report
        .first_annihilating_multiple()
        .expect("Hasse's theorem should produce an annihilating multiple");

    assert!(report.interval().contains(annihilating_multiple));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(annihilating_multiple).expect("small-field Hasse multiple should fit in u64"),
    ));
}

#[test]
fn finite_group_surfaces_match_the_short_companion_on_a_reducible_curve() {
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(curve.order(), conversion.target().order());
    assert_eq!(curve.exponent(), conversion.target().exponent());
    assert_eq!(
        curve.order_distribution(),
        conversion.target().order_distribution()
    );
    assert_eq!(
        curve.group_structure(),
        conversion.target().group_structure()
    );
    assert_eq!(
        curve.describe_group_structure(),
        conversion.target().describe_group_structure()
    );
    assert_eq!(
        curve
            .frobenius_trace()
            .expect("Montgomery curve should yield a Frobenius trace"),
        conversion
            .target()
            .frobenius_trace()
            .expect("short companion should yield a Frobenius trace"),
    );
}

#[test]
fn finite_field_wrappers_remain_honest_about_short_reduction_limits() {
    let curve = f3_curve();

    assert!(matches!(
        curve.group_order_by(FiniteFieldGroupOrderStrategy::Auto),
        Err(CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(3u8)
    ));
    assert!(matches!(
        curve.frobenius_trace_by(FiniteFieldGroupOrderStrategy::Schoof),
        Err(CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(3u8)
    ));
}

#[test]
fn group_order_by_small_field_exhaustive_is_native_in_characteristic_three() {
    let curve = f3_curve();

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Exhaustive)
        .expect("native exhaustive small-field route should succeed");

    assert!(matches!(report, GroupOrderReport::ExhaustiveTrace(_)));
    assert_eq!(report.route(), GroupOrderRoute::Exhaustive);
    assert_eq!(report.curve_order(), BigUint::from(4u8));
}

#[test]
fn group_order_by_small_field_auto_falls_back_to_exhaustive_when_short_routes_are_unavailable() {
    let curve = f3_curve();

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Auto)
        .expect("auto route should fall back to exhaustive");

    assert!(matches!(report, GroupOrderReport::ExhaustiveTrace(_)));
    assert_eq!(report.route(), GroupOrderRoute::Exhaustive);
}

#[test]
fn group_order_by_small_field_quadratic_character_still_matches_the_short_companion_when_supported()
{
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(
        curve
            .group_order_by_small_field(SmallFieldGroupOrderStrategy::QuadraticCharacter)
            .expect("Montgomery wrapper should succeed")
            .curve_order(),
        conversion
            .target()
            .group_order_by_small_field(SmallFieldGroupOrderStrategy::QuadraticCharacter)
            .expect("short companion route should succeed")
            .curve_order(),
    );
}

#[test]
fn point_order_by_exhaustive_is_native_in_characteristic_three() {
    let curve = f3_curve();
    let point = curve
        .point(F3::from_i64(2), F3::one())
        .expect("sample point should lie on the curve");

    let report = curve
        .point_order_by(&point, PointOrderStrategy::Exhaustive)
        .expect("native exhaustive point-order route should succeed");

    let PointOrderReport::Exhaustive(report) = report else {
        panic!("expected the exhaustive point-order variant");
    };
    assert_eq!(report.exact_order(), &num_bigint::BigUint::from(4u8));
}

#[test]
fn point_order_by_hasse_interval_naive_uses_the_requested_small_field_route() {
    let curve = f3_curve();
    let point = curve
        .point(F3::from_i64(2), F3::one())
        .expect("sample point should lie on the curve");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy: SmallFieldGroupOrderStrategy::Exhaustive,
            },
        )
        .expect("native Hasse-interval route should succeed");

    let PointOrderReport::HasseIntervalNaive(report) = report else {
        panic!("expected the HasseIntervalNaive variant");
    };
    assert_eq!(report.exact_order(), &num_bigint::BigUint::from(4u8));
    assert_eq!(
        report.group_order_report().route(),
        GroupOrderRoute::Exhaustive
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_montgomery_finite_group_surfaces_are_self_consistent_in_characteristic_three(
        curve in arb_nonsingular_montgomery_curve::<crate::fields::Fp3>(CurveStrategyConfig::default()),
    ) {
        let structure = curve.group_structure();
        let trace = curve.frobenius_trace().expect("small finite curve should yield a trace");

        prop_assert_eq!(curve.check_group_axioms(), Ok(()));
        prop_assert_eq!(curve.order(), curve.points().len());
        prop_assert_eq!(structure.order, curve.order());
        prop_assert_eq!(structure.exponent, curve.exponent());
        prop_assert_eq!(trace.curve_order(), BigUint::from(curve.order() as u64));
        prop_assert!(curve.verify_hasse_bound().expect("Hasse bound check should succeed").holds());
    }

    #[test]
    fn property_montgomery_finite_group_surfaces_match_the_short_companion_in_characteristic_greater_than_three(
        (curve, point) in arb_montgomery_curve_and_point::<crate::fields::Fp5>(CurveStrategyConfig::default()),
    ) {
        let conversion = curve
            .conversion_to_short_weierstrass()
            .expect("characteristic five should support the short reduction");
        let short_point = conversion
            .map_source_point(&point)
            .expect("sampled point should transport to short");
        let general_point_order = curve.point_order(&point).expect("sampled point should have an order");
        let short_point_order = conversion.target().point_order(&short_point).expect("transported point should have an order");
        let hasse_interval = HasseInterval::for_field::<F5>().expect("valid Hasse interval");

        prop_assert_eq!(curve.order(), conversion.target().order());
        prop_assert_eq!(curve.exponent(), conversion.target().exponent());
        prop_assert_eq!(curve.order_distribution(), conversion.target().order_distribution());
        prop_assert_eq!(curve.group_structure(), conversion.target().group_structure());
        prop_assert_eq!(curve.describe_group_structure(), conversion.target().describe_group_structure());
        prop_assert_eq!(general_point_order, short_point_order);
        prop_assert_eq!(curve.point_has_exact_order(&point, general_point_order), Ok(true));
        prop_assert_eq!(
            curve
                .frobenius_trace()
                .expect("Montgomery curve should yield a Frobenius trace"),
            conversion
                .target()
                .frobenius_trace()
                .expect("short companion should yield a Frobenius trace"),
        );
        prop_assert_eq!(
            curve
                .find_annihilating_multiple_in_hasse_interval_naive(&point)
                .expect("naive Hasse search should succeed")
                .interval()
                .clone(),
            hasse_interval,
        );
    }
}
