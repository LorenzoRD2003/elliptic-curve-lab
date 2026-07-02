use super::shared::{F5, F13, f5_curve, f13_denominator_curve};
use crate::elliptic_curves::{
    CurveError,
    traits::{AffineCurveModel, GroupCurveModel, HasProjectiveModel, ProjectiveGroupCurveModel},
    twisted_edwards::projective::{
        ExtendedTwistedEdwardsPoint, TwistedEdwardsProjectiveOperationCost,
        TwistedEdwardsProjectiveOperationKind,
    },
};
use crate::fields::FieldError;
use crate::fields::traits::*;

#[test]
fn projective_negation_matches_affine_negation_after_recovery() {
    let curve = f5_curve();
    let affine = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");
    let projective = curve
        .to_projective(&affine)
        .expect("projective lift should succeed");

    let recovered = curve
        .to_affine_projective(&curve.neg_projective(&projective))
        .expect("projective negation should recover affinely");

    assert_eq!(recovered, curve.neg(&affine));
}

#[test]
fn projective_addition_matches_affine_addition_on_sample_points() {
    let curve = f5_curve();
    let left = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should be on the curve");
    let right = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");

    let projective_sum = curve
        .add_projective(
            &curve.to_projective(&left).expect("lift should succeed"),
            &curve.to_projective(&right).expect("lift should succeed"),
        )
        .expect("projective addition should succeed");
    let recovered = curve
        .to_affine_projective(&projective_sum)
        .expect("projective sum should recover affinely");

    assert_eq!(
        recovered,
        curve
            .add(&left, &right)
            .expect("affine addition should succeed")
    );
}

#[test]
fn projective_doubling_matches_affine_doubling() {
    let curve = f5_curve();
    let affine = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should be on the curve");
    let projective = curve
        .to_projective(&affine)
        .expect("projective lift should succeed");

    let recovered = curve
        .to_affine_projective(
            &curve
                .double_projective(&projective)
                .expect("projective doubling should succeed"),
        )
        .expect("projective double should recover affinely");

    assert_eq!(
        recovered,
        curve
            .double(&affine)
            .expect("affine doubling should succeed")
    );
}

#[test]
fn projective_scalar_multiplication_matches_affine_scalar_multiplication() {
    let curve = f5_curve();
    let affine = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");
    let projective = curve
        .to_projective(&affine)
        .expect("projective lift should succeed");

    let recovered = curve
        .to_affine_projective(
            &curve
                .mul_scalar_projective(&projective, 3)
                .expect("projective scalar multiplication should succeed"),
        )
        .expect("projective scalar multiple should recover affinely");

    assert_eq!(
        recovered,
        curve
            .mul_scalar(&affine, 3)
            .expect("affine scalar multiplication should succeed")
    );
}

#[test]
fn mixed_projective_addition_matches_affine_addition_on_sample_points() {
    let curve = f5_curve();
    let left = curve
        .point(F5::one(), F5::zero())
        .expect("sample point should be on the curve");
    let right = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");

    let mixed_sum = curve
        .mixed_add_projective(
            &curve.to_projective(&left).expect("lift should succeed"),
            &right,
        )
        .expect("mixed projective addition should succeed");
    let recovered = curve
        .to_affine_projective(&mixed_sum)
        .expect("mixed projective sum should recover affinely");

    assert_eq!(
        recovered,
        curve
            .add(&left, &right)
            .expect("affine addition should succeed")
    );
}

#[test]
fn projective_group_operations_reject_points_outside_the_curve() {
    let curve = f5_curve();
    let valid = curve
        .to_projective(
            &curve
                .point(F5::one(), F5::zero())
                .expect("sample point should be on the curve"),
        )
        .expect("lift should succeed");
    let invalid = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(2),
        F5::from_i64(2),
        F5::one(),
        F5::zero(),
    );

    assert_eq!(
        curve.add_projective(&valid, &invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        curve.double_projective(&invalid),
        Err(CurveError::PointNotOnCurve)
    );
    assert_eq!(
        curve.mul_scalar_projective(&invalid, 3),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn projective_cost_reports_expose_native_extended_operation_counts() {
    let add =
        TwistedEdwardsProjectiveOperationCost::for_kind(TwistedEdwardsProjectiveOperationKind::Add);
    let mixed = TwistedEdwardsProjectiveOperationCost::for_kind(
        TwistedEdwardsProjectiveOperationKind::MixedAdd,
    );
    let scalar = TwistedEdwardsProjectiveOperationCost::for_scalar_mul(13);

    assert_eq!(add.affine_additions(), 0);
    assert_eq!(add.affine_doublings(), 0);
    assert_eq!(add.representation_cost().multiplications(), 9);
    assert_eq!(mixed.representation_cost().multiplications(), 8);
    assert_eq!(scalar.representation_cost().multiplications(), 45);
    assert_eq!(scalar.representation_cost().squarings(), 9);
}

#[test]
fn native_projective_doubling_can_leave_the_affine_chart_honestly() {
    let curve = f13_denominator_curve();
    let affine = curve
        .point(F13::from_i64(4), F13::from_i64(3))
        .expect("sample point should be on the curve");
    let projective = curve
        .to_projective(&affine)
        .expect("projective lift should succeed");
    let doubled = curve
        .double_projective(&projective)
        .expect("native projective doubling should stay defined");

    assert!(curve.is_projective_point_on_curve(&doubled));
    assert!(F13::is_zero(doubled.z()));
    assert_eq!(
        curve.to_affine_projective(&doubled),
        Err(CurveError::Field(FieldError::DivisionByZero))
    );
}

#[test]
fn affine_wrappers_still_surface_honest_denominator_failures() {
    let curve = f13_denominator_curve();
    let affine = curve
        .point(F13::from_i64(4), F13::from_i64(3))
        .expect("sample point should be on the curve");

    assert_eq!(
        curve.double(&affine),
        Err(CurveError::Field(FieldError::DivisionByZero))
    );
}
