use super::shared::{F3, F5, F7, f3_curve, f5_curve, f7_curve};
use crate::elliptic_curves::{
    AffinePoint, CurveError,
    traits::{
        AffineCurveModel, CurveModel, CurveModelConversion, FiniteGroupCurveModel, GroupCurveModel,
    },
};
use crate::fields::traits::Field;

#[test]
fn montgomery_negation_uses_the_native_affine_involution() {
    let curve = f5_curve();
    let point = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the curve");

    assert_eq!(
        curve.neg(&point),
        curve
            .point(F5::from_i64(2), F5::from_i64(3))
            .expect("the inverse should lie on the curve")
    );
    assert_eq!(curve.neg(&curve.identity()), curve.identity());
}

#[test]
fn montgomery_addition_matches_the_short_companion_transport() {
    let curve = f5_curve();
    let left = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("left point should lie on the curve");
    let right = curve
        .point(F5::from_i64(3), F5::from_i64(2))
        .expect("right point should lie on the curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let expected = conversion
        .map_target_point(
            &conversion
                .target()
                .add(
                    &conversion
                        .map_source_point(&left)
                        .expect("left point should transport to short"),
                    &conversion
                        .map_source_point(&right)
                        .expect("right point should transport to short"),
                )
                .expect("short companion addition should succeed"),
        )
        .expect("short sum should transport back");

    assert_eq!(curve.add(&left, &right), Ok(expected));
}

#[test]
fn montgomery_doubling_matches_the_short_companion_transport() {
    let curve = f5_curve();
    let point = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should lie on the curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let expected = conversion
        .map_target_point(
            &conversion
                .target()
                .double(
                    &conversion
                        .map_source_point(&point)
                        .expect("point should transport to short"),
                )
                .expect("short companion doubling should succeed"),
        )
        .expect("short double should transport back");

    assert_eq!(curve.double(&point), Ok(expected));
}

#[test]
fn montgomery_small_torsion_helpers_work_with_the_native_affine_group_law() {
    let curve = f5_curve();
    let point = curve
        .point(F5::zero(), F5::zero())
        .expect("sample point should lie on the curve");

    assert_eq!(curve.double(&point), Ok(curve.identity()));
    assert_eq!(curve.point_has_exact_order(&point, 2), Ok(true));
    assert!(curve.is_torsion_point(&point, 2));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(curve.identity()));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn montgomery_group_law_handles_characteristic_three_natively() {
    let curve = f3_curve();
    let point = curve
        .point(F3::from_i64(2), F3::one())
        .expect("sample point should lie on the curve");
    let inverse = curve
        .point(F3::from_i64(2), F3::from_i64(2))
        .expect("the inverse point should lie on the curve");
    let doubled = curve
        .point(F3::zero(), F3::zero())
        .expect("the doubled point should lie on the curve");

    assert_eq!(curve.neg(&point), inverse);
    assert_eq!(curve.add(&point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.double(&point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(doubled));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn montgomery_group_operations_reject_points_outside_the_curve() {
    let curve = f7_curve();
    let valid = curve
        .point(F7::from_i64(2), F7::one())
        .expect("point should lie on the curve");
    let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

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
    assert_eq!(
        curve.mul_scalar_signed(&invalid, -3),
        Err(CurveError::PointNotOnCurve)
    );
}
