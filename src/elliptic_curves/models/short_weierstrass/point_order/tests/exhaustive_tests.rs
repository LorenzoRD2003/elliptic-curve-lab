use crate::elliptic_curves::{
    AffinePoint, CurveError,
    short_weierstrass::point_order::{
        ExhaustivePointOrderReport, PointOrderReport, PointOrderStrategy, PointOrderStrategyKind,
    },
    traits::AffineCurveModel,
};
use crate::fields::traits::Field;

use super::shared::{F7, bu, f7_curve};

#[test]
fn exhaustive_route_reports_the_exact_order() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("valid point");

    let report = curve
        .point_order_by(&point, PointOrderStrategy::Exhaustive)
        .expect("exhaustive route should succeed");

    assert_eq!(report.strategy_kind(), PointOrderStrategyKind::Exhaustive);
    assert_eq!(report.exact_order(), &bu(6));
    assert_eq!(
        report,
        PointOrderReport::Exhaustive(ExhaustivePointOrderReport::new(bu(6)))
    );
}

#[test]
fn point_order_by_rejects_points_outside_the_curve() {
    let curve = f7_curve();
    let invalid = AffinePoint::<F7>::new(F7::from_i64(1), F7::from_i64(1));

    assert_eq!(
        curve.point_order_by(&invalid, PointOrderStrategy::Exhaustive),
        Err(CurveError::PointNotOnCurve)
    );
}
