use super::shared::{bu, f7_curve, f7_point, point_index, sampler_from_indices};
use crate::elliptic_curves::{
    frobenius::group_order::{GroupOrderRoute, SmallFieldGroupOrderStrategy},
    short_weierstrass::group_exponent::{GroupExponentReport, GroupExponentStrategy},
    short_weierstrass::point_order::{PointOrderReport, PointOrderStrategy},
};

#[test]
fn report_strategy_preserves_the_selected_route() {
    let curve = f7_curve();
    let mut sampler = |_| Some(0usize);

    let exhaustive_report = curve
        .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
        .expect("exhaustive exponent route should succeed");

    assert_eq!(
        exhaustive_report.strategy(),
        GroupExponentStrategy::Exhaustive
    );

    let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(2, 1))]);
    let random_report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("random-point route should succeed");

    assert_eq!(
        random_report.strategy(),
        GroupExponentStrategy::RandomPoints {
            max_samples: 1,
            point_order_strategy: PointOrderStrategy::Exhaustive,
        }
    );
}

#[test]
fn group_exponent_by_exhaustive_recovers_the_exact_small_group_exponent() {
    let curve = f7_curve();
    let mut sampler = |_| Some(0usize);

    let report = curve
        .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
        .expect("exhaustive exponent route should succeed");

    assert_eq!(report.strategy(), GroupExponentStrategy::Exhaustive);
    assert_eq!(report.exponent_lower_bound(), &bu(6));
    assert_eq!(report.exact_exponent(), Some(&bu(6)));

    let GroupExponentReport::Exhaustive(exponent) = report else {
        panic!("expected the exhaustive route to preserve its variant");
    };
    assert_eq!(exponent, bu(6));
}

#[test]
fn group_exponent_by_random_points_accumulates_lcms_of_sampled_point_orders() {
    let curve = f7_curve();
    let point_order_two = f7_point(6, 0);
    let point_order_six = f7_point(2, 1);
    let mut sampler = sampler_from_indices(vec![
        point_index(&curve, &point_order_two),
        point_index(&curve, &point_order_six),
    ]);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 2,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("random-point accumulation should succeed");

    assert_eq!(
        report.strategy(),
        GroupExponentStrategy::RandomPoints {
            max_samples: 2,
            point_order_strategy: PointOrderStrategy::Exhaustive,
        }
    );
    assert_eq!(report.exponent_lower_bound(), &bu(6));
    assert_eq!(report.exact_exponent(), None);

    let GroupExponentReport::RandomPoints(report) = report else {
        panic!("expected the random-point route to preserve its variant");
    };

    assert_eq!(report.samples_requested(), 2);
    assert_eq!(report.samples_taken(), 2);
    assert!(report.completed_requested_samples());
    assert_eq!(
        report.point_order_strategy(),
        &PointOrderStrategy::Exhaustive
    );
    assert_eq!(
        report.strategy(),
        GroupExponentStrategy::RandomPoints {
            max_samples: 2,
            point_order_strategy: PointOrderStrategy::Exhaustive,
        }
    );
    assert_eq!(report.exponent_lower_bound(), &bu(6));
    assert_eq!(report.steps().len(), 2);
    assert_eq!(report.steps()[0].point(), &point_order_two);
    assert_eq!(report.steps()[0].accumulated_lcm(), &bu(2));
    assert_eq!(report.steps()[1].point(), &point_order_six);
    assert_eq!(report.steps()[1].accumulated_lcm(), &bu(6));
}

#[test]
fn group_exponent_by_random_points_preserves_the_point_order_route() {
    let curve = f7_curve();
    let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(2, 1))]);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::HasseIntervalNaive {
                    group_order_strategy: SmallFieldGroupOrderStrategy::Auto,
                },
            },
            &mut sampler,
        )
        .expect("Hasse-driven point-order route should compose into exponent accumulation");

    let GroupExponentReport::RandomPoints(report) = report else {
        panic!("expected the random-point route to preserve its variant");
    };

    match report.steps()[0].point_order_report() {
        PointOrderReport::HasseIntervalNaive(step_report) => {
            assert_eq!(
                step_report.group_order_report().route(),
                GroupOrderRoute::QuadraticCharacter
            );
            assert_eq!(step_report.exact_order(), &bu(6));
        }
        other => panic!("expected HasseIntervalNaive point-order report, got {other:?}"),
    }
}

#[test]
fn group_exponent_by_random_points_reports_early_sampler_exhaustion_honestly() {
    let curve = f7_curve();
    let mut sampler = sampler_from_indices(vec![point_index(&curve, &f7_point(6, 0))]);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 3,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("random-point accumulation should succeed");

    let GroupExponentReport::RandomPoints(report) = report else {
        panic!("expected the random-point route to preserve its variant");
    };

    assert_eq!(report.samples_requested(), 3);
    assert_eq!(report.samples_taken(), 1);
    assert!(!report.completed_requested_samples());
    assert_eq!(report.exponent_lower_bound(), &bu(2));
}
