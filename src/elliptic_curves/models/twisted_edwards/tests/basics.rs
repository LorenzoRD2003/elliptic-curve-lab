use super::shared::{F2, F5, f5_curve};
use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant},
};
use crate::fields::traits::Field;

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
    let curve = f5_curve();

    assert!(F5::eq(curve.a(), &F5::one()));
    assert!(F5::eq(curve.d(), &F5::from_i64(2)));
    assert_eq!(
        curve.to_equation_string(),
        "(1 (mod 5))x^2 + y^2 = 1 + (2 (mod 5))x^2y^2"
    );
}

#[test]
fn invariants_match_known_small_example() {
    let curve = f5_curve();

    assert!(F5::eq(&curve.c4(), &F5::from_i64(3)));
    assert!(F5::eq(&curve.c6(), &F5::one()));
    assert!(F5::eq(&curve.discriminant(), &F5::from_i64(2)));
    assert!(F5::eq(&curve.j_invariant(), &F5::one()));
}

#[test]
fn invariants_satisfy_weierstrass_identity() {
    let curve = f5_curve();

    let left = F5::sub(&F5::cube(&curve.c4()), &F5::square(&curve.c6()));
    let right = F5::mul(&F5::from_i64(1728), &curve.discriminant());

    assert!(F5::eq(&left, &right));
}

#[test]
fn same_curve_has_same_j_invariant() {
    let left = f5_curve();
    let right = f5_curve();

    assert!(left.has_same_j_invariant(&right));
}

#[test]
fn curve_model_identity_is_the_finite_point_zero_one() {
    let curve = f5_curve();
    let identity = curve.identity();

    assert!(curve.is_identity(&identity));
    assert!(curve.contains(&identity));
}

#[test]
fn affine_infinity_is_not_on_the_twisted_edwards_model() {
    let curve = f5_curve();

    assert!(!curve.contains(&AffinePoint::<F5>::infinity()));
}

#[test]
fn point_constructor_accepts_points_on_the_curve_and_rejects_points_off_it() {
    let curve = f5_curve();

    assert!(curve.point(F5::zero(), F5::one()).is_ok());
    assert_eq!(
        curve.point(F5::one(), F5::one()),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_method() {
    let curve = f5_curve();

    assert!(F5::eq(
        &<TwistedEdwardsCurve<F5> as HasJInvariant>::j_invariant(&curve),
        &curve.j_invariant()
    ));
}
