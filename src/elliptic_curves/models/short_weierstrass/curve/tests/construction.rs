use crate::elliptic_curves::traits::{AffineCurveModel, CurveModel};
use crate::elliptic_curves::{CurveError, ShortWeierstrassCurve};
use crate::fields::Q;
use crate::fields::traits::*;
use num_bigint::BigUint;

use super::shared::{F2, F3, F5, F7, f7_curve, q};

#[test]
fn constructor_rejects_characteristics_two_and_three() {
    assert!(matches!(
        ShortWeierstrassCurve::<F2>::new(F2::zero(), F2::one()),
        Err(CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(2u8),
    ));
    assert!(matches!(
        ShortWeierstrassCurve::<F3>::new(F3::zero(), F3::one()),
        Err(CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(3u8),
    ));
}

#[test]
fn constructor_rejects_singular_coefficients() {
    assert!(matches!(
        ShortWeierstrassCurve::<F5>::new(F5::zero(), F5::zero()),
        Err(CurveError::SingularCurve),
    ));
}

#[test]
fn characteristic_zero_fields_are_allowed() {
    let curve = ShortWeierstrassCurve::<Q>::new(q(-1, 1), q(0, 1)).expect("non-singular curve");
    let point = curve
        .point(q(0, 1), q(0, 1))
        .expect("point should lie on the curve");

    assert!(curve.contains(&point));
    assert!(Q::has_characteristic(0));
}

#[test]
fn point_constructor_accepts_valid_affine_coordinates() {
    let curve = f7_curve();

    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");

    assert!(matches!(
        point,
        crate::elliptic_curves::AffinePoint::Finite { .. }
    ));
}

#[test]
fn point_constructor_rejects_invalid_affine_coordinates() {
    let curve =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve");

    assert!(matches!(
        curve.point(F7::from_i64(2), F7::from_i64(2)),
        Err(CurveError::PointNotOnCurve)
    ));
}
