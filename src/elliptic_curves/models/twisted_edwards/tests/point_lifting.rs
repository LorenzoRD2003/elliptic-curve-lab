use super::shared::{F5, F13, f5_curve, f13_denominator_curve};
use crate::elliptic_curves::{
    AffinePoint,
    traits::{CurveModel, LiftXCoordinate, LiftedPoints},
};
use crate::fields::traits::Field;

#[test]
fn lift_x_at_zero_returns_the_identity_and_its_inverse() {
    let curve = f5_curve();

    match curve.lift_x(F5::zero()).expect("lifting should succeed") {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.is_identity(&left) || curve.is_identity(&right));
            assert!(
                left == AffinePoint::new(F5::zero(), F5::one())
                    || right == AffinePoint::new(F5::zero(), F5::one())
            );
            assert!(
                left == AffinePoint::new(F5::zero(), F5::from_i64(-1))
                    || right == AffinePoint::new(F5::zero(), F5::from_i64(-1))
            );
        }
        other => panic!("expected two lifted points above x = 0, got {other:?}"),
    }
}

#[test]
fn lift_x_returns_no_point_when_the_fiber_denominator_vanishes() {
    let curve = f13_denominator_curve();

    assert_eq!(
        curve.lift_x(F13::one()).expect("lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn point_from_x_reports_none_when_no_affine_point_exists() {
    let curve = f13_denominator_curve();

    assert_eq!(curve.point_from_x(F13::one()), Ok(None));
}
