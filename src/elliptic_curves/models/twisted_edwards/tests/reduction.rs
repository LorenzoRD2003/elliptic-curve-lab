use super::shared::{F5, F13, f5_curve, f13_curve};
use crate::elliptic_curves::{GeneralWeierstrassCurve, MontgomeryCurve, ShortWeierstrassCurve};
use crate::fields::traits::Field;

#[test]
fn twisted_edwards_to_montgomery_matches_known_small_example() {
    let curve = f5_curve();
    let montgomery = curve.as_montgomery();

    assert!(F5::eq(montgomery.a(), &F5::from_i64(4)));
    assert!(F5::eq(montgomery.b(), &F5::one()));
}

#[test]
fn twisted_edwards_to_montgomery_roundtrip_preserves_coefficients() {
    let twisted = f13_curve();
    let roundtrip = twisted.as_montgomery().as_twisted_edwards();

    assert!(F13::eq(twisted.a(), roundtrip.a()));
    assert!(F13::eq(twisted.d(), roundtrip.d()));
}

#[test]
fn montgomery_to_twisted_edwards_roundtrip_preserves_coefficients() {
    let montgomery = MontgomeryCurve::<F13>::new(F13::from_i64(3), F13::from_i64(2))
        .expect("sample Montgomery curve should be non-singular");
    let roundtrip = montgomery.as_twisted_edwards().as_montgomery();

    assert!(F13::eq(montgomery.a(), roundtrip.a()));
    assert!(F13::eq(montgomery.b(), roundtrip.b()));
}

#[test]
fn direct_and_composed_montgomery_j_invariants_agree() {
    let twisted = f13_curve();
    let montgomery = twisted.as_montgomery();

    assert!(F13::eq(&twisted.j_invariant(), &montgomery.j_invariant()));
}

#[test]
fn composed_general_weierstrass_conversion_preserves_j_invariant() {
    let twisted = f13_curve();
    let general: GeneralWeierstrassCurve<F13> = twisted.as_montgomery().as_general_weierstrass();

    assert!(F13::eq(&twisted.j_invariant(), &general.j_invariant()));
}

#[test]
fn composed_short_weierstrass_conversion_preserves_j_invariant() {
    let twisted = f5_curve();
    let short: ShortWeierstrassCurve<F5> = twisted
        .as_montgomery()
        .try_as_short_weierstrass()
        .expect("characteristic 5 should support the Montgomery short companion");

    assert!(F5::eq(&twisted.j_invariant(), &short.j_invariant()));
}
