use proptest::prelude::*;
use std::collections::HashSet;

use crate::elliptic_curves::frobenius::HasseInterval;
use crate::elliptic_curves::frobenius::hasse::search::HasseIntervalSearchCurveModel;
use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModel, CurveModelConversion, EnumerableCurveModel,
    FiniteGroupCurveModel, FrobeniusTraceCurveModel, GroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve};
use crate::fields::{Fp, traits::Field};
use crate::isogenies::kernel::IsogenyKernel;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::{
    arb_general_weierstrass_curve_and_point, arb_nonsingular_general_weierstrass_curve,
};

type F2 = Fp<2>;
type F3 = Fp<3>;
type F5 = Fp<5>;

#[test]
fn enumerable_curve_model_lists_the_expected_points_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    let finite_points = curve.finite_points();
    let all_points = curve.points();

    assert_eq!(finite_points, vec![AffinePoint::new(F2::one(), F2::zero())]);
    assert_eq!(all_points.len(), 2);
    assert_eq!(all_points[0], curve.identity());
}

#[test]
fn finite_group_curve_model_reports_the_expected_structure_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let structure = curve.group_structure();

    assert_eq!(curve.order(), 2);
    assert_eq!(curve.exponent(), 2);
    assert_eq!(curve.order_distribution().get(&1), Some(&1));
    assert_eq!(curve.order_distribution().get(&2), Some(&1));
    assert_eq!(
        structure,
        crate::elliptic_curves::traits::FiniteAbelianGroupStructure {
            order: 2,
            exponent: 2,
            cyclic: true,
            invariant_factors: None,
        }
    );
}

#[test]
fn frobenius_trace_curve_model_recovers_the_expected_trace_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let trace = curve
        .frobenius_trace()
        .expect("small finite curve should yield a Frobenius trace");

    assert_eq!(trace.curve_order(), 2);
    assert_eq!(trace.field_order(), 2);
    assert_eq!(trace.trace(), 1);
    assert!(
        curve
            .verify_hasse_bound()
            .expect("Hasse bound verification should succeed")
            .holds()
    );
}

#[test]
fn hasse_search_finds_an_annihilating_multiple_for_a_general_curve_point() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::from_i64(2), F5::one())
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
fn isogeny_kernel_cyclic_works_for_general_weierstrass_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let generator = curve
        .point(F5::from_i64(2), F5::one())
        .expect("sample point should lie on the curve");
    let kernel =
        IsogenyKernel::cyclic(&curve, &generator).expect("point should generate a cyclic subgroup");

    let expected = vec![curve.identity(), generator.clone()];

    assert_eq!(curve.point_order(&generator), Some(2));
    assert_eq!(kernel.points(), expected.as_slice());
    assert_eq!(kernel.order(), 2);
}

#[test]
fn explicit_isogeny_kernel_validation_works_for_general_weierstrass_subgroups() {
    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let generator = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let cyclic_kernel =
        IsogenyKernel::cyclic(&curve, &generator).expect("generator should define a cyclic kernel");
    let explicit_kernel = IsogenyKernel::new(
        &curve,
        cyclic_kernel
            .points()
            .iter()
            .cloned()
            .collect::<HashSet<_>>(),
    )
    .expect("the explicitly re-listed cyclic subgroup should stay valid");

    assert_eq!(explicit_kernel.order(), cyclic_kernel.order());
    for point in cyclic_kernel.points() {
        assert!(explicit_kernel.contains(point));
    }
}

#[test]
fn finite_group_surfaces_match_the_short_companion_on_a_reducible_curve() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
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
            .expect("general curve should yield a Frobenius trace"),
        conversion
            .target()
            .frobenius_trace()
            .expect("short companion should yield a Frobenius trace"),
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_general_weierstrass_finite_group_surfaces_are_self_consistent_in_characteristic_two(
        curve in arb_nonsingular_general_weierstrass_curve::<2>(CurveStrategyConfig::default()),
    ) {
        let structure = curve.group_structure();
        let trace = curve.frobenius_trace().expect("small finite curve should yield a trace");

        prop_assert_eq!(curve.check_group_axioms(), Ok(()));
        prop_assert_eq!(curve.order(), curve.points().len());
        prop_assert_eq!(structure.order, curve.order());
        prop_assert_eq!(structure.exponent, curve.exponent());
        prop_assert_eq!(trace.curve_order(), curve.order() as u64);
        prop_assert!(curve.verify_hasse_bound().expect("Hasse bound check should succeed").holds());
    }

    #[test]
    fn property_general_weierstrass_finite_group_surfaces_are_self_consistent_in_characteristic_three(
        curve in arb_nonsingular_general_weierstrass_curve::<3>(CurveStrategyConfig::default()),
    ) {
        let structure = curve.group_structure();
        let trace = curve.frobenius_trace().expect("small finite curve should yield a trace");

        prop_assert_eq!(curve.check_group_axioms(), Ok(()));
        prop_assert_eq!(curve.order(), curve.points().len());
        prop_assert_eq!(structure.order, curve.order());
        prop_assert_eq!(structure.exponent, curve.exponent());
        prop_assert_eq!(trace.curve_order(), curve.order() as u64);
        prop_assert!(curve.verify_hasse_bound().expect("Hasse bound check should succeed").holds());
    }

    #[test]
    fn property_general_weierstrass_finite_group_surfaces_match_the_short_companion_in_characteristic_greater_than_three(
        (curve, point) in arb_general_weierstrass_curve_and_point::<5>(CurveStrategyConfig::default()),
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
                .expect("general curve should yield a Frobenius trace"),
            conversion
                .target()
                .frobenius_trace()
                .expect("short companion should yield a Frobenius trace"),
        );

        let general_report = curve
            .find_annihilating_multiple_in_interval_naive(&point, hasse_interval.clone())
            .expect("general curve should search the Hasse interval");
        let short_report = conversion
            .target()
            .find_annihilating_multiple_in_interval_naive(&short_point, hasse_interval)
            .expect("short companion should search the Hasse interval");

        prop_assert_eq!(
            general_report.first_annihilating_multiple(),
            short_report.first_annihilating_multiple()
        );

        let general_kernel = IsogenyKernel::cyclic(&curve, &point)
            .expect("sampled point should generate a cyclic subgroup");
        let short_kernel = IsogenyKernel::cyclic(conversion.target(), &short_point)
            .expect("transported point should generate a cyclic subgroup");

        prop_assert_eq!(general_kernel.order(), short_kernel.order());
    }
}
