use crate::fields::traits::*;
use crypto_bigint::{U64, U256, const_prime_monty_params};

use crate::elliptic_curves::{
    GeneralWeierstrassCurve, MontgomeryCurve, ShortWeierstrassCurve, TwistedEdwardsCurve,
    traits::{
        AffineCurveModel, CurveModel, GroupCurveModel, HasProjectiveModel, LiftXCoordinate,
        LiftedPoints, ProjectiveGroupCurveModel,
    },
};
use crate::fields::Fp;

const_prime_monty_params!(P17Prime, U64, "0000000000000011", 3);
const_prime_monty_params!(
    P256Prime,
    U256,
    "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
    6
);

type F17 = Fp<P17Prime, { U64::LIMBS }>;
type Fp256 = Fp<P256Prime, { U256::LIMBS }>;

fn f17(value: i64) -> <F17 as crate::fields::traits::Field>::Elem {
    F17::from_i64(value)
}

fn fp256(value: i64) -> <Fp256 as Field>::Elem {
    Fp256::from_i64(value)
}

fn fp256_div(numerator: i64, denominator: i64) -> <Fp256 as Field>::Elem {
    Fp256::div(&fp256(numerator), &fp256(denominator)).expect("denominator should be non-zero")
}

#[test]
fn short_weierstrass_supports_montgomery_prime_fields() {
    let curve = ShortWeierstrassCurve::<F17>::new(f17(2), f17(2))
        .expect("non-singular short-Weierstrass curve over F17");
    let point = curve.point(f17(5), f17(1)).expect("known F17 point");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(f17(5)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));

    let doubled = curve.double(&point).expect("doubling should work");
    assert_eq!(curve.add(&point, &point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 3), curve.add(&doubled, &point));

    let projective = curve
        .to_projective(&point)
        .expect("affine lift should work");
    let projective_double = curve
        .double_projective(&projective)
        .expect("projective doubling should work");
    assert_eq!(curve.to_affine_projective(&projective_double), Ok(doubled));
}

#[test]
fn short_weierstrass_accepts_a_large_montgomery_prime() {
    let curve = ShortWeierstrassCurve::<Fp256>::new(fp256(2), fp256(3))
        .expect("non-singular short-Weierstrass curve over P-256");
    let point = curve
        .point(fp256(3), fp256(6))
        .expect("(3, 6) lies on y^2 = x^3 + 2x + 3");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(fp256(3)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));
}

#[test]
fn general_weierstrass_supports_montgomery_prime_fields() {
    let curve = GeneralWeierstrassCurve::<F17>::new(f17(0), f17(0), f17(0), f17(2), f17(2))
        .expect("short-form general-Weierstrass curve over F17");
    let point = curve.point(f17(5), f17(1)).expect("known F17 point");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(f17(5)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));

    let doubled = curve.double(&point).expect("doubling should work");
    assert_eq!(curve.add(&point, &point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 3), curve.add(&doubled, &point));

    let projective = curve
        .to_projective(&point)
        .expect("affine lift should work");
    let projective_double = curve
        .double_projective(&projective)
        .expect("projective doubling should work");
    assert_eq!(curve.to_affine_projective(&projective_double), Ok(doubled));
}

#[test]
fn general_weierstrass_accepts_a_large_montgomery_prime() {
    let curve =
        GeneralWeierstrassCurve::<Fp256>::new(fp256(0), fp256(0), fp256(0), fp256(2), fp256(3))
            .expect("short-form general-Weierstrass curve over P-256");
    let point = curve
        .point(fp256(3), fp256(6))
        .expect("(3, 6) lies on the large-prime curve");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(fp256(3)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));
}

#[test]
fn montgomery_curves_support_montgomery_prime_fields() {
    let curve = MontgomeryCurve::<F17>::new(f17(0), f17(1))
        .expect("non-singular Montgomery curve over F17");
    let point = curve
        .point(f17(1), f17(6))
        .expect("(1, 6) lies on y^2 = x^3 + x over F17");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(f17(1)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));

    let doubled = curve.double(&point).expect("doubling should work");
    assert_eq!(curve.add(&point, &point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 3), curve.add(&doubled, &point));
}

#[test]
fn montgomery_curves_accept_a_large_montgomery_prime() {
    let curve = MontgomeryCurve::<Fp256>::new(fp256(5), fp256_div(75, 16))
        .expect("non-singular Montgomery curve over P-256");
    let point = curve
        .point(fp256(3), fp256(4))
        .expect("(3, 4) was used to choose B");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(fp256(3)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));
}

#[test]
fn twisted_edwards_supports_montgomery_prime_fields() {
    let curve = TwistedEdwardsCurve::<F17>::new(f17(1), f17(6))
        .expect("non-singular twisted-Edwards curve over F17");
    let point = curve
        .point(f17(2), f17(3))
        .expect("(2, 3) lies on x^2 + y^2 = 1 + 6x^2y^2 over F17");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(f17(2)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));

    let doubled = curve.double(&point).expect("doubling should work");
    assert_eq!(curve.add(&point, &point), Ok(doubled.clone()));
    assert_eq!(curve.mul_scalar(&point, 3), curve.add(&doubled, &point));

    let projective = curve
        .to_projective(&point)
        .expect("affine lift should work");
    let projective_double = curve
        .double_projective(&projective)
        .expect("projective doubling should work");
    assert_eq!(curve.to_affine_projective(&projective_double), Ok(doubled));
}

#[test]
fn twisted_edwards_accepts_a_large_montgomery_prime() {
    let curve = TwistedEdwardsCurve::<Fp256>::new(fp256(1), fp256_div(1, 3))
        .expect("non-singular twisted-Edwards curve over P-256");
    let point = curve
        .point(fp256(2), fp256(3))
        .expect("(2, 3) was used to choose d = 1/3");

    assert!(curve.contains(&point));
    assert!(matches!(
        curve.lift_x(fp256(2)),
        Ok(LiftedPoints::TwoPoints(_, _))
    ));
}
