use crate::elliptic_curves::traits::{CurveModel, LiftXCoordinate};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};
use crate::fields::{Q, traits::Field};

use super::shared::{F7, f7_curve, q};

#[test]
fn point_from_x_returns_one_point_when_rhs_has_a_square_root() {
    let curve = f7_curve();

    let point = curve
        .point_from_x(F7::from_i64(2))
        .expect("x = 2 should lift to a point");

    super::shared::assert_contains(&curve, &point);
    assert!(matches!(point, AffinePoint::Finite { .. }));
}

#[test]
fn point_from_x_returns_none_when_rhs_is_not_a_square() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(curve.point_from_x(F7::from_i64(0)).is_none());
}

#[test]
fn points_from_x_returns_both_points_when_they_are_distinct() {
    let curve = f7_curve();

    let (left, right) = curve
        .points_from_x(F7::from_i64(2))
        .expect("x = 2 should lift to two points");

    super::shared::assert_contains(&curve, &left);
    super::shared::assert_contains(&curve, &right);
    assert_ne!(left, right);
}

#[test]
fn points_from_x_repeats_the_point_when_the_square_root_is_zero() {
    let curve = f7_curve();

    let (left, right) = curve
        .points_from_x(F7::from_i64(6))
        .expect("x = 6 should give y = 0");

    assert_eq!(left, right);
    super::shared::assert_contains(&curve, &left);
}

#[test]
fn points_from_x_works_over_q_when_an_exact_root_exists() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

    let (left, right) = curve
        .points_from_x(q(1, 1))
        .expect("x = 1 should give y = 0 in Q");

    assert_eq!(left, right);
    assert!(curve.contains(&left));
}
