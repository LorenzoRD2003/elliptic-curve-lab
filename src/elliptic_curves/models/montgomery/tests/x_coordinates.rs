use super::shared::F5;
use crate::elliptic_curves::MontgomeryXzPoint;
use crate::fields::traits::Field;

#[test]
fn xz_points_compare_equal_up_to_projective_rescaling() {
    let left = MontgomeryXzPoint::<F5>::new(F5::from_i64(2), F5::one());
    let right = MontgomeryXzPoint::<F5>::new(F5::from_i64(4), F5::from_i64(2));

    assert_eq!(left, right);
    assert!(!left.has_same_representative_as(&right));
}

#[test]
fn xz_point_normalization_rescales_finite_points_to_z_equal_one() {
    let point = MontgomeryXzPoint::<F5>::new(F5::from_i64(4), F5::from_i64(2));

    let normalized = point
        .normalize()
        .expect("nonzero Z should admit affine x-recovery");

    assert_eq!(
        normalized,
        MontgomeryXzPoint::<F5>::from_affine_x(F5::from_i64(2))
    );
    assert!(normalized.is_normalized());
}

#[test]
fn xz_points_roundtrip_affine_x_when_z_is_invertible() {
    let point = MontgomeryXzPoint::<F5>::new(F5::from_i64(4), F5::from_i64(2));

    assert_eq!(
        point
            .to_affine_x()
            .expect("nonzero Z should recover one affine x"),
        Some(F5::from_i64(2))
    );
}

#[test]
fn xz_infinity_is_explicit_and_already_normalized() {
    let infinity = MontgomeryXzPoint::<F5>::infinity();

    assert!(infinity.is_infinity());
    assert!(infinity.is_normalized());
    assert_eq!(
        infinity
            .to_affine_x()
            .expect("infinity should recover cleanly"),
        None
    );
}

#[test]
fn finite_xz_points_with_zero_z_fail_affine_recovery_and_normalization() {
    let point = MontgomeryXzPoint::<F5>::new(F5::one(), F5::zero());

    assert!(point.to_affine_x().is_err());
    assert!(point.normalize().is_err());
}
