use crate::elliptic_curves::{
    GeneralWeierstrassCurve,
    traits::{CurveModel, LiftXCoordinate, LiftedPoints},
};
use crate::fields::{ComplexApprox, Q, traits::Field};

use super::shared::{F2, F4, F5, c};

#[test]
fn lift_x_over_an_odd_characteristic_curve_returns_two_points_when_the_fiber_is_split() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let lifted = curve
        .lift_x(F5::zero())
        .expect("finite-field odd-characteristic lifting should succeed");

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
fn lift_x_over_an_odd_characteristic_curve_returns_no_point_when_the_fiber_is_empty() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    assert_eq!(
        curve
            .lift_x(F5::from_i64(3))
            .expect("finite-field odd-characteristic lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_one_point_when_b_is_zero() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    assert_eq!(
        curve
            .lift_x(F2::one())
            .expect("characteristic-two lifting should succeed"),
        LiftedPoints::OnePoint(crate::elliptic_curves::AffinePoint::<F2>::new(
            F2::one(),
            F2::zero()
        ))
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_no_point_when_the_artin_schreier_equation_is_unsolvable() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    assert_eq!(
        curve
            .lift_x(F2::zero())
            .expect("characteristic-two lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_in_characteristic_two_returns_two_points_when_the_artin_schreier_equation_is_solvable() {
    let curve = GeneralWeierstrassCurve::<F4>::new(
        F4::zero(),
        F4::zero(),
        F4::one(),
        F4::zero(),
        F4::zero(),
    )
    .expect("non-singular curve in characteristic two");

    let lifted = curve
        .lift_x(F4::one())
        .expect("characteristic-two lifting should succeed");

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
fn lift_x_over_q_returns_two_rational_points_when_the_fiber_is_split() {
    let curve =
        GeneralWeierstrassCurve::<Q>::new(Q::one(), Q::one(), Q::one(), Q::one(), Q::zero())
            .expect("non-singular curve over Q");

    let lifted = curve
        .lift_x(Q::zero())
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
fn lift_x_over_q_returns_no_point_when_the_quadratic_has_no_rational_root() {
    let curve =
        GeneralWeierstrassCurve::<Q>::new(Q::one(), Q::one(), Q::one(), Q::one(), Q::zero())
            .expect("non-singular curve over Q");

    assert_eq!(
        curve
            .lift_x(Q::from_i64(3))
            .expect("rational lifting should succeed"),
        LiftedPoints::NoPoint
    );
}

#[test]
fn lift_x_over_complex_approx_returns_two_points() {
    let curve = GeneralWeierstrassCurve::<ComplexApprox>::new(
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(1.0, 0.0),
        c(0.0, 0.0),
    )
    .expect("non-singular curve over ComplexApprox");

    let lifted = curve
        .lift_x(c(3.0, 0.0))
        .expect("complex lifting should succeed");

    match lifted {
        LiftedPoints::TwoPoints(left, right) => {
            assert!(curve.contains(&left));
            assert!(curve.contains(&right));
            assert_ne!(left, right);
        }
        other => panic!("expected two lifted points over ComplexApprox, got {other:?}"),
    }
}
