use super::shared::{
    alternate_f7_curve, alternate_f7_point, bu, f7_curve, f7_point, point_index,
    sampler_from_indices,
};
use crate::elliptic_curves::{
    CurveError,
    frobenius::group_order::GroupOrderStrategy,
    short_weierstrass::group_exponent::{
        ExponentLowerBoundGroupOrderVerification, GroupExponentReport, GroupExponentStrategy,
    },
    short_weierstrass::point_order::PointOrderStrategy,
};

#[test]
fn verify_exponent_lower_bound_by_group_order_can_report_non_unique_hasse_multiples() {
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
        .expect("random-point route should succeed");
    let GroupExponentReport::RandomPoints(accumulation) = report else {
        panic!("expected random-point accumulation");
    };

    let verification = curve
        .verify_exponent_lower_bound_by_group_order(&accumulation, GroupOrderStrategy::Auto)
        .expect("verification should succeed");

    assert_eq!(verification.exponent_lower_bound(), &bu(6));
    assert_eq!(verification.verified_group_order(), None);
}

#[test]
fn verify_exponent_lower_bound_by_group_order_rejects_reports_from_other_curves() {
    let curve = f7_curve();
    let other_curve = alternate_f7_curve();
    let mut sampler =
        sampler_from_indices(vec![point_index(&other_curve, &alternate_f7_point(2, 2))]);

    let report = other_curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("other-curve accumulation should succeed");
    let GroupExponentReport::RandomPoints(accumulation) = report else {
        panic!("expected random-point accumulation");
    };

    let err = curve
        .verify_exponent_lower_bound_by_group_order(&accumulation, GroupOrderStrategy::Auto)
        .expect_err("cross-curve accumulation should be rejected");
    assert_eq!(err, CurveError::PointNotOnCurve);
}

#[test]
fn verification_report_exposes_its_inputs_and_unique_multiple() {
    let curve = f7_curve();
    let verification = curve.verify_exponent_lower_bound_by_group_order_report(
        bu(6),
        curve
            .group_order_by(GroupOrderStrategy::Auto)
            .expect("group order should succeed"),
    );

    let cloned: ExponentLowerBoundGroupOrderVerification = verification.clone();
    assert_eq!(cloned.exponent_lower_bound(), &bu(6));
    assert_eq!(cloned.group_order_report().curve_order(), 6);
    assert_eq!(cloned.verified_group_order(), None);
}
