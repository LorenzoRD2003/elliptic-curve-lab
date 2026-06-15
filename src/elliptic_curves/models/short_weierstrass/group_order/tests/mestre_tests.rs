use proptest::prelude::*;

use super::shared::{
    f241_curve, genuine_twist_curve, max_order_point_index, sampler_covering_each_curve_by_index,
};
use crate::elliptic_curves::{
    frobenius::group_order::{GroupOrderReport, GroupOrderStrategy, MestreConfig},
    traits::EnumerableCurveModel,
};
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};

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
        mestre_report.resolved_side_label(),
        "original curve" | "quadratic twist"
    ));
    assert!(matches!(
        mestre_report.resolved_side_group_order_candidate(),
        candidate if candidate == mestre_report.curve_order()
            || candidate == mestre_report.twist_curve_order()
    ));
    assert!(mestre_report.step_count() > 0);
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
            report.resolved_side_label(),
            "original curve" | "quadratic twist"
        ));
        prop_assert!(matches!(
            report.resolved_side_group_order_candidate(),
            candidate if candidate == report.curve_order()
                || candidate == report.twist_curve_order()
        ));
        prop_assert!(report.step_count() > 0);
        prop_assert!(report.step_count() <= max_iterations);
    }
}
