use crate::elliptic_curves::traits::{CurveModel, LiftXCoordinate, LiftedPoints};
use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve};
use crate::fields::Q;

use super::shared::{F7, f7_curve, q};

#[test]
fn lift_x_returns_two_points_when_the_square_roots_are_distinct() {
    let curve = f7_curve();

    let lifted = curve
        .lift_x(F7::from_i64(2))
        .expect("lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            super::shared::assert_contains(&curve, &left);
            super::shared::assert_contains(&curve, &right);
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points, got {other:?}"),
    }
}

#[test]
fn lift_x_returns_no_point_when_rhs_is_not_a_square() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert_eq!(
        curve
            .lift_x(F7::from_i64(0))
            .expect("lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn point_from_x_returns_one_point_when_the_fiber_is_nonempty() {
    let curve = f7_curve();

    let point = curve
        .point_from_x(F7::from_i64(2))
        .expect("lifting should succeed")
        .expect("x = 2 should lift to a point");

    super::shared::assert_contains(&curve, &point);
    assert!(matches!(point, AffinePoint::Finite { .. }));
}

#[test]
fn points_from_x_repeats_the_point_when_the_fiber_has_one_point() {
    let curve = f7_curve();

    let (left, right) = curve
        .points_from_x(F7::from_i64(6))
        .expect("lifting should succeed")
        .expect("x = 6 should give y = 0");

    super::shared::assert_contains(&curve, &left);
    super::shared::assert_contains(&curve, &right);
    assert_eq!(left, right);
}

#[test]
fn points_from_x_works_over_q_when_an_exact_root_exists() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");

    let (left, right) = curve
        .points_from_x(q(1, 1))
        .expect("lifting should succeed")
        .expect("x = 1 should give y = 0 in Q");

    assert_eq!(left, right);
    assert!(curve.contains(&left));
}
