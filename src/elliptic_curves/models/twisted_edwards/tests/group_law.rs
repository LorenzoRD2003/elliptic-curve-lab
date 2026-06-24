use crate::elliptic_curves::{
    AffinePoint, CurveError,
    traits::{AffineCurveModel, CurveModel, GroupCurveModel},
};
use crate::fields::{FieldError, traits::Field};

use super::shared::{F5, F13, f5_curve, f13_denominator_curve};

#[test]
fn twisted_edwards_negation_uses_the_native_affine_involution() {
    let curve = f5_curve();
    let point = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");

    assert_eq!(
        curve.neg(&point),
        curve
            .point(F5::from_i64(-1), F5::zero())
            .expect("the inverse should lie on the curve")
    );
    assert_eq!(curve.neg(&curve.identity()), curve.identity());
}

#[test]
fn twisted_edwards_addition_matches_montgomery_on_the_common_open_subset() {
    let curve = f5_curve();
    let left = curve
        .point(F5::one(), F5::zero())
        .expect("left point should lie on the twisted-Edwards curve");
    let right = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("right point should lie on the twisted-Edwards curve");
    let montgomery = curve.as_montgomery();

    let expected = montgomery
        .try_point_to_twisted_edwards_open(
            &montgomery
                .add(
                    &curve
                        .try_point_to_montgomery_open(&left)
                        .expect("left point should transport to Montgomery"),
                    &curve
                        .try_point_to_montgomery_open(&right)
                        .expect("right point should transport to Montgomery"),
                )
                .expect("Montgomery addition should succeed on the common open"),
        )
        .expect("the Montgomery sum should stay in the common open");

    assert_eq!(curve.add(&left, &right), Ok(expected));
}

#[test]
fn twisted_edwards_addition_handles_identity_and_inverse_cases() {
    let curve = f5_curve();
    let point = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");
    let inverse = curve.neg(&point);

    assert_eq!(curve.add(&curve.identity(), &point), Ok(point.clone()));
    assert_eq!(curve.add(&point, &curve.identity()), Ok(point.clone()));
    assert_eq!(curve.add(&point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.add(&inverse, &point), Ok(curve.identity()));
}

#[test]
fn twisted_edwards_doubling_and_scalar_multiplication_work_on_a_small_example() {
    let curve = f5_curve();
    let point = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");

    assert_eq!(curve.double(&point), curve.add(&point, &point));
    assert_eq!(
        curve.mul_scalar(&point, 3),
        curve.add(
            &curve.double(&point).expect("doubling should succeed"),
            &point
        )
    );
}

#[test]
fn twisted_edwards_group_operations_reject_points_outside_the_curve() {
    let curve = f5_curve();
    let valid = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should lie on the twisted-Edwards curve");
    let invalid = AffinePoint::<F5>::new(F5::one(), F5::one());

    assert_eq!(
        curve.add(&valid, &invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(curve.double(&invalid), Err(CurveError::PointNotOnCurve));
    assert_eq!(
        curve.sub(&valid, &invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        curve.mul_scalar(&invalid, 3),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn twisted_edwards_generic_affine_formulas_report_zero_denominators_honestly() {
    let curve = f13_denominator_curve();
    let left = curve
        .point(F13::from_i64(4), F13::from_i64(3))
        .expect("sample point should lie on the twisted-Edwards curve");

    assert_eq!(
        curve.double(&left),
        Err(CurveError::Field(FieldError::DivisionByZero))
    );
}
