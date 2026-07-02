use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    traits::{CurveModelConversion, CurveModelConversionError},
};
use crate::fields::traits::*;
use num_bigint::BigUint;

use super::shared::{F5, F7, f3_curve, f5_curve, f7_scaled_curve};

#[test]
fn reduction_rejects_characteristic_three() {
    let curve = f3_curve();

    assert!(matches!(
        curve.conversion_to_short_weierstrass(),
        Err(CurveModelConversionError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(3u8),
    ));
}

#[test]
fn reduction_produces_expected_short_companion_for_a_small_example() {
    let curve = f5_curve();

    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    assert_eq!(conversion.source(), &curve);
    assert!(F5::eq(conversion.target().a(), &F5::from_i64(4)));
    assert!(F5::eq(conversion.target().b(), &F5::from_i64(4)));
}

#[test]
fn reduction_produces_expected_short_companion_when_b_is_not_one() {
    let curve = f7_scaled_curve();

    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic seven should support the reduction");

    assert!(F7::eq(conversion.target().a(), &F7::from_i64(3)));
    assert!(F7::eq(conversion.target().b(), &F7::one()));
}

#[test]
fn try_as_short_weierstrass_matches_the_reduction_companion() {
    let curve = f5_curve();

    let from_helper = curve
        .try_as_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let from_reduction = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction")
        .target()
        .clone();

    assert_eq!(from_helper, from_reduction);
}

#[test]
fn try_from_montgomery_reference_matches_the_reduction_companion() {
    let curve = f5_curve();

    let from_try_from = ShortWeierstrassCurve::try_from(&curve)
        .expect("characteristic five should support the reduction");
    let from_reduction = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction")
        .target()
        .clone();

    assert_eq!(from_try_from, from_reduction);
}

#[test]
fn transporting_infinity_between_montgomery_and_short_models_is_stable() {
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");

    let short_infinity = conversion
        .map_source_point(&AffinePoint::<F5>::Infinity)
        .expect("infinity should transport to the short companion");
    let montgomery_infinity = conversion
        .map_target_point(&AffinePoint::<F5>::Infinity)
        .expect("infinity should transport back to the Montgomery model");

    assert_eq!(short_infinity, AffinePoint::<F5>::Infinity);
    assert_eq!(montgomery_infinity, AffinePoint::<F5>::Infinity);
}

#[test]
fn transporting_a_montgomery_point_to_short_matches_the_expected_coordinates() {
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let montgomery_point = AffinePoint::<F5>::new(F5::zero(), F5::zero());

    let image = conversion
        .map_source_point(&montgomery_point)
        .expect("point should transport to the short companion");

    assert_eq!(image, AffinePoint::<F5>::new(F5::from_i64(2), F5::zero()));
}

#[test]
fn transporting_short_and_montgomery_points_roundtrips() {
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the reduction");
    let montgomery_point = AffinePoint::<F5>::new(F5::from_i64(2), F5::from_i64(2));
    let short_point = conversion
        .map_source_point(&montgomery_point)
        .expect("point should transport to the short companion");

    let montgomery_roundtrip = conversion
        .map_target_point(&short_point)
        .expect("short point should transport back to the Montgomery model");
    let short_roundtrip = conversion
        .map_source_point(&montgomery_roundtrip)
        .expect("transporting back again should still succeed");

    assert_eq!(montgomery_roundtrip, montgomery_point);
    assert_eq!(short_roundtrip, short_point);
}

#[test]
fn conversion_reports_invalid_source_and_target_points_honestly() {
    let curve = f5_curve();
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the conversion");
    let bad_montgomery_point = AffinePoint::<F5>::new(F5::zero(), F5::one());
    let bad_short_point = AffinePoint::<F5>::new(F5::one(), F5::one());

    assert_eq!(
        conversion.map_source_point(&bad_montgomery_point),
        Err(CurveModelConversionError::PointNotOnSource)
    );
    assert_eq!(
        conversion.map_target_point(&bad_short_point),
        Err(CurveModelConversionError::PointNotOnTarget)
    );
}
