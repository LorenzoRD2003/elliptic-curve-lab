use crate::elliptic_curves::{
    frobenius::group_order::{GroupOrderRoute, SmallFieldGroupOrderStrategy},
    short_weierstrass::point_order::{
        HasseIntervalPointOrderReport, PointOrderReport, PointOrderStrategy, PointOrderStrategyKind,
    },
    traits::AffineCurveModel,
};
use crate::numerics::NormalizedPrimePowerFactorization;

use super::shared::{F7, bu, f7_curve};

#[test]
fn known_multiple_route_preserves_the_prime_peeling_report() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("valid point");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::FromKnownMultiple {
                multiple: bu(6),
                factorization: vec![(bu(2), 1), (bu(3), 1)],
            },
        )
        .expect("known-multiple route should succeed");

    assert_eq!(
        report.strategy_kind(),
        PointOrderStrategyKind::FromKnownMultiple
    );
    assert_eq!(report.exact_order(), &bu(6));
}

#[test]
fn hasse_interval_route_records_group_order_search_and_prime_peeling() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("valid point");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy: SmallFieldGroupOrderStrategy::Auto,
            },
        )
        .expect("Hasse-interval route should succeed");

    assert_eq!(
        report.strategy_kind(),
        PointOrderStrategyKind::HasseIntervalNaive
    );

    let PointOrderReport::HasseIntervalNaive(report) = report else {
        panic!("expected Hasse-interval route to preserve its variant");
    };

    let HasseIntervalPointOrderReport {
        group_order_report,
        multiple_search,
        order_from_multiple,
    } = &*report;

    assert_eq!(
        group_order_report.route(),
        GroupOrderRoute::QuadraticCharacter
    );
    assert_eq!(multiple_search.first_annihilating_multiple(), Some(&bu(6)));
    assert_eq!(order_from_multiple.exact_order(), &bu(6));
}

#[test]
fn hasse_interval_route_reuses_trusted_factorization_from_found_multiple() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(6), F7::from_i64(0))
        .expect("valid point");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy: SmallFieldGroupOrderStrategy::Auto,
            },
        )
        .expect("Hasse-interval route should succeed");

    let PointOrderReport::HasseIntervalNaive(report) = report else {
        panic!("expected Hasse-interval route to preserve its variant");
    };

    let found_multiple = report
        .multiple_search()
        .first_annihilating_multiple()
        .expect("the Hasse route should have found one annihilating multiple")
        .clone();
    let expected_factorization = NormalizedPrimePowerFactorization::factor(&found_multiple)
        .expect("the found multiple should factor")
        .into_factors();
    let from_report: Vec<_> = report
        .order_from_multiple()
        .steps()
        .iter()
        .map(|step| (step.prime().clone(), step.exponent_in_multiple()))
        .collect();

    assert_eq!(from_report, expected_factorization);
    assert_eq!(report.exact_order(), &bu(2));
}
