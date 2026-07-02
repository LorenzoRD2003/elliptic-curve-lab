use super::shared::{F3, F5, F7, f3_curve, f5_curve, f7_curve, q};
use crate::elliptic_curves::{
    AffinePoint, MontgomeryCurve,
    traits::{CurveModel, LiftXCoordinate, LiftedPoints},
};
use crate::fields::Q;
use crate::fields::traits::*;

#[test]
fn lift_x_over_a_small_prime_field_returns_two_points_when_the_fiber_is_split() {
    let curve = f5_curve();

    let lifted = curve
        .lift_x(F5::from_i64(2))
        .expect("finite-field lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points, got {other:?}"),
    }
}

#[test]
fn lift_x_over_a_small_prime_field_returns_one_point_when_y_is_zero() {
    let curve = f5_curve();

    assert_eq!(
        curve
            .lift_x(F5::zero())
            .expect("finite-field lifting should succeed"),
        LiftedPoints::OnePoint(AffinePoint::<F5>::new(F5::zero(), F5::zero()))
    );
}

#[test]
fn lift_x_over_a_small_prime_field_returns_no_point_when_the_fiber_is_empty() {
    let curve = f5_curve();

    assert_eq!(
        curve
            .lift_x(F5::one())
            .expect("finite-field lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_in_characteristic_three_still_works() {
    let curve = f3_curve();

    let lifted = curve
        .lift_x(F3::from_i64(2))
        .expect("characteristic-three lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points in characteristic three, got {other:?}"),
    }
}

#[test]
fn lift_x_over_q_returns_two_rational_points_when_the_fiber_is_split() {
    let curve = MontgomeryCurve::<Q>::new(q(7, 1), q(1, 1)).expect("non-singular curve over Q");

    let lifted = curve
        .lift_x(q(1, 1))
        .expect("rational lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points over Q, got {other:?}"),
    }
}

#[test]
fn lift_x_over_q_returns_no_point_when_the_fiber_has_no_rational_root() {
    let curve = MontgomeryCurve::<Q>::new(q(7, 1), q(1, 1)).expect("non-singular curve over Q");

    assert_eq!(
        curve
            .lift_x(q(2, 1))
            .expect("rational lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn point_from_x_and_points_from_x_reuse_the_lifted_fiber() {
    let curve = f7_curve();

    let one_point = curve
        .point_from_x(F7::from_i64(2))
        .expect("lifting helpers should succeed")
        .expect("x = 2 should lift");
    let pair = curve
        .points_from_x(F7::from_i64(2))
        .expect("lifting helpers should succeed")
        .expect("x = 2 should lift");

    assert!(curve.contains(&one_point));
    assert!(curve.contains(&pair.0));
    assert!(curve.contains(&pair.1));
}
