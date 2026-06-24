use crate::elliptic_curves::{CurveError, TwistedEdwardsCurve};
use crate::fields::{Fp, traits::Field};

type F2 = Fp<2>;
type F5 = Fp<5>;

#[test]
fn constructor_rejects_characteristic_two() {
    assert!(matches!(
        TwistedEdwardsCurve::<F2>::new(F2::one(), F2::zero()),
        Err(CurveError::UnsupportedCharacteristic { characteristic: 2 })
    ));
}

#[test]
fn constructor_rejects_zero_a_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::zero(), F5::one()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn constructor_rejects_zero_d_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::one(), F5::zero()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn constructor_rejects_equal_a_and_d_as_singular() {
    assert!(matches!(
        TwistedEdwardsCurve::<F5>::new(F5::one(), F5::one()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn valid_curve_exposes_coefficients_and_equation_string() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(F5::eq(curve.a(), &F5::one()));
    assert!(F5::eq(curve.d(), &F5::from_i64(2)));
    assert_eq!(
        curve.to_equation_string(),
        "(1 (mod 5))x^2 + y^2 = 1 + (2 (mod 5))x^2y^2"
    );
}

#[test]
fn invariants_match_known_small_example() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(F5::eq(&curve.c4(), &F5::from_i64(3)));
    assert!(F5::eq(&curve.c6(), &F5::one()));
    assert!(F5::eq(&curve.discriminant(), &F5::from_i64(2)));
    assert!(F5::eq(&curve.j_invariant(), &F5::one()));
}

#[test]
fn invariants_satisfy_weierstrass_identity() {
    let curve = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    let left = F5::sub(&F5::cube(&curve.c4()), &F5::square(&curve.c6()));
    let right = F5::mul(&F5::from_i64(1728), &curve.discriminant());

    assert!(F5::eq(&left, &right));
}

#[test]
fn same_curve_has_same_j_invariant() {
    let left = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");
    let right = TwistedEdwardsCurve::<F5>::new(F5::one(), F5::from_i64(2))
        .expect("sample twisted-Edwards curve should be non-singular");

    assert!(left.has_same_j_invariant(&right));
}
