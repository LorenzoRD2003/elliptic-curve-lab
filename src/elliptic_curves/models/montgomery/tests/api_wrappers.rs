use super::shared::{F3, F5, F7, f5_curve, f7_scaled_curve};
use crate::elliptic_curves::{
    GeneralWeierstrassCurve, MontgomeryCurve, ShortWeierstrassCurve,
    traits::CurveModelConversionError,
};
use crate::fields::traits::*;
use num_bigint::BigUint;

#[test]
fn as_general_weierstrass_matches_the_expected_direct_embedding() {
    let curve = f7_scaled_curve();

    let general = curve.as_general_weierstrass();

    assert!(F7::eq(general.a1(), &F7::zero()));
    assert!(F7::eq(general.a2(), &F7::from_i64(5)));
    assert!(F7::eq(general.a3(), &F7::zero()));
    assert!(F7::eq(general.a4(), &F7::from_i64(2)));
    assert!(F7::eq(general.a6(), &F7::zero()));
}

#[test]
fn from_montgomery_reference_matches_as_general_weierstrass() {
    let curve = f5_curve();

    assert_eq!(
        GeneralWeierstrassCurve::from(&curve),
        curve.as_general_weierstrass()
    );
}

#[test]
fn short_to_montgomery_roundtrips_back_to_the_same_short_curve() {
    let short_curve =
        ShortWeierstrassCurve::<F5>::new(F5::from_i64(4), F5::from_i64(4)).expect("valid curve");

    let montgomery = MontgomeryCurve::try_from(&short_curve)
        .expect("this small short curve should admit a Montgomery model");
    let short_roundtrip = montgomery
        .try_as_short_weierstrass()
        .expect("the reconstructed Montgomery curve should reduce back to short");

    assert_eq!(short_roundtrip, short_curve);
}

#[test]
fn general_to_montgomery_roundtrips_on_a_small_reducible_curve() {
    let montgomery = f5_curve();
    let general = montgomery.as_general_weierstrass();

    let recovered = MontgomeryCurve::try_from(&general)
        .expect("the embedded general curve should recover a Montgomery model");
    let recovered_short = recovered
        .try_as_short_weierstrass()
        .expect("the recovered Montgomery model should reduce to short");
    let original_short = montgomery
        .try_as_short_weierstrass()
        .expect("the original Montgomery model should reduce to short");

    assert_eq!(recovered_short, original_short);
}

#[test]
fn try_as_montgomery_matches_try_from_on_general_curves() {
    let general = f5_curve().as_general_weierstrass();

    assert_eq!(
        general
            .try_as_montgomery()
            .expect("embedded general curve should admit a Montgomery model"),
        MontgomeryCurve::try_from(&general)
            .expect("embedded general curve should admit a Montgomery model"),
    );
}

#[test]
fn general_to_montgomery_fails_honestly_in_characteristic_three() {
    let general = GeneralWeierstrassCurve::<F3>::new(
        F3::zero(),
        F3::zero(),
        F3::zero(),
        F3::one(),
        F3::zero(),
    )
    .expect("non-singular characteristic-three curve");

    assert!(matches!(
        MontgomeryCurve::try_from(&general),
        Err(CurveModelConversionError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(3u8)
    ));
}
