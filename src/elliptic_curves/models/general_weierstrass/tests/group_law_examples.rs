use crate::elliptic_curves::{
    GeneralWeierstrassCurve,
    traits::{
        AffineCurveModel, CurveModel, CurveModelConversion, FiniteGroupCurveModel, GroupCurveModel,
    },
};
use crate::fields::{Fp, traits::Field};

use super::shared::{F2, F5};

#[test]
fn general_weierstrass_negation_uses_the_model_specific_involution() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::zero(), F5::from_i64(4))
        .expect("sample point should lie on the curve");

    assert_eq!(
        curve.neg(&point),
        curve
            .point(F5::zero(), F5::zero())
            .expect("the inverse should lie on the curve")
    );
    assert_ne!(curve.neg(&point), point.neg());
}

#[test]
fn general_weierstrass_addition_matches_the_short_companion_transport() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let left = curve
        .point(F5::zero(), F5::zero())
        .expect("left point should lie on the curve");
    let right = curve
        .point(F5::from_i64(2), F5::one())
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
fn general_weierstrass_doubling_matches_the_short_companion_transport() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::zero(), F5::zero())
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
fn general_weierstrass_small_torsion_helpers_work_with_the_native_affine_group_law() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let point = curve
        .point(F5::from_i64(2), F5::one())
        .expect("sample point should lie on the curve");

    assert_eq!(curve.double(&point), Ok(curve.identity()));
    assert_eq!(curve.point_has_exact_order(&point, 2), Ok(true));
    assert!(curve.is_torsion_point(&point, 2));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(curve.identity()));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn general_weierstrass_group_law_handles_characteristic_two_natively() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let point = curve
        .point(F2::one(), F2::zero())
        .expect("sample point should lie on the curve");

    assert_eq!(curve.add(&point, &point), Ok(curve.identity()));
    assert_eq!(curve.neg(&point), point);
    assert_eq!(curve.double(&point), Ok(curve.identity()));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(curve.identity()));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}

#[test]
fn general_weierstrass_group_law_handles_characteristic_three_natively() {
    type F3 = Fp<3>;

    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let point = curve
        .point(F3::zero(), F3::zero())
        .expect("sample point should lie on the curve");
    let inverse = curve
        .point(F3::zero(), F3::from_i64(2))
        .expect("the inverse point should lie on the curve");
    let doubled = curve
        .point(F3::from_i64(2), F3::one())
        .expect("the doubled point should lie on the curve");

    assert_eq!(curve.neg(&point), inverse);
    assert_eq!(curve.add(&point, &inverse), Ok(curve.identity()));
    assert_eq!(curve.double(&point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 2), Ok(doubled));
    assert_eq!(curve.check_group_axioms(), Ok(()));
}
