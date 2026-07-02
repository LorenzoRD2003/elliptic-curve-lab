use super::shared::{F2, F3, F5, F7, q};
use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve, MontgomeryCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant},
};
use crate::fields::Q;
use crate::fields::traits::*;
use num_bigint::BigUint;

#[test]
fn constructor_rejects_characteristic_two() {
    assert!(matches!(
        MontgomeryCurve::<F2>::new(F2::zero(), F2::one()),
        Err(CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(2u8)
    ));
}

#[test]
fn constructor_allows_characteristic_three_when_the_curve_is_nonsingular() {
    let curve = MontgomeryCurve::<F3>::new(F3::zero(), F3::one())
        .expect("Montgomery model should be valid in characteristic three");

    assert!(F3::eq(curve.a(), &F3::zero()));
    assert!(F3::eq(curve.b(), &F3::one()));
    assert!(F3::eq(&curve.discriminant(), &F3::from_i64(2)));
}

#[test]
fn constructor_rejects_zero_b_and_singular_a_squared_minus_four() {
    assert!(matches!(
        MontgomeryCurve::<F5>::new(F5::zero(), F5::zero()),
        Err(CurveError::SingularCurve)
    ));
    assert!(matches!(
        MontgomeryCurve::<F5>::new(F5::from_i64(2), F5::one()),
        Err(CurveError::SingularCurve)
    ));
}

#[test]
fn accessors_and_invariants_match_a_small_prime_field_example() {
    let curve =
        MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::one()).expect("valid Montgomery curve");

    assert!(F7::eq(curve.a(), &F7::from_i64(3)));
    assert!(F7::eq(curve.b(), &F7::one()));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c4(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(4)));
    assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(2)));
    assert!(F7::eq(&curve.rhs_value(&F7::zero()), &F7::zero()));
}

#[test]
fn montgomery_invariants_satisfy_the_classical_relation() {
    let curve =
        MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::one()).expect("valid Montgomery curve");

    let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
    let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

    assert!(F7::eq(&left, &right));
}

#[test]
fn invariants_match_the_equivalent_general_weierstrass_model() {
    let montgomery =
        MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid curve");
    let general = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::div(montgomery.a(), montgomery.b()).expect("B != 0"),
        F7::zero(),
        F7::div(&F7::one(), &F7::square(montgomery.b())).expect("B != 0"),
        F7::zero(),
    )
    .expect("the scaled Montgomery equation should define a valid general model");

    assert!(F7::eq(&montgomery.c4(), &general.c4()));
    assert!(F7::eq(&montgomery.c6(), &general.c6()));
    assert!(F7::eq(&montgomery.discriminant(), &general.discriminant()));
    assert!(F7::eq(&montgomery.j_invariant(), &general.j_invariant()));
}

#[test]
fn discriminant_uses_the_expected_b_to_the_sixth_scaling() {
    let curve = MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid curve");

    assert!(F7::eq(&curve.c4(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(4)));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
}

#[test]
fn j_invariant_matches_a_classical_exact_example_over_q() {
    let curve = MontgomeryCurve::<Q>::new(q(0, 1), q(1, 1)).expect("non-singular curve");

    assert!(Q::eq(&curve.discriminant(), &q(-64, 1)));
    assert!(Q::eq(&curve.c4(), &q(-48, 1)));
    assert!(Q::eq(&curve.c6(), &q(0, 1)));
    assert!(Q::eq(&curve.j_invariant(), &q(1728, 1)));
}

#[test]
fn display_and_debug_surface_the_montgomery_equation() {
    let curve = MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid curve");

    let display = curve.to_string();

    assert!(display.contains("y^2 = x^3 + "));
    assert!(display.contains("x^2 + x"));
    assert!(format!("{curve:?}").contains("MontgomeryCurve"));
}

#[test]
fn curve_model_identity_and_membership_helpers_work_for_montgomery() {
    let curve = MontgomeryCurve::<F5>::new(F5::zero(), F5::one()).expect("valid Montgomery curve");
    let finite_point = AffinePoint::<F5>::new(F5::zero(), F5::zero());
    let infinity = AffinePoint::<F5>::infinity();
    let off_curve_point = AffinePoint::<F5>::new(F5::zero(), F5::one());

    assert_eq!(curve.identity(), infinity);
    assert!(curve.contains(&finite_point));
    assert!(curve.is_on_curve_nonzero(&finite_point));
    assert!(curve.contains(&infinity));
    assert!(curve.is_identity(&infinity));
    assert!(!curve.is_on_curve_nonzero(&infinity));
    assert!(!curve.contains(&off_curve_point));
    assert!(!curve.is_on_curve_nonzero(&off_curve_point));
}

#[test]
fn affine_curve_model_point_accepts_valid_montgomery_points() {
    let curve = MontgomeryCurve::<F5>::new(F5::zero(), F5::one()).expect("valid Montgomery curve");

    let point = curve
        .point(F5::zero(), F5::zero())
        .expect("the affine point should lie on the curve");

    assert_eq!(point, AffinePoint::<F5>::new(F5::zero(), F5::zero()));
    assert!(curve.contains(&point));
}

#[test]
fn affine_curve_model_point_rejects_invalid_montgomery_points() {
    let curve = MontgomeryCurve::<F5>::new(F5::zero(), F5::one()).expect("valid Montgomery curve");

    assert_eq!(
        curve.point(F5::zero(), F5::one()),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_method() {
    let curve =
        MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::one()).expect("valid Montgomery curve");

    assert!(F7::eq(
        &HasJInvariant::j_invariant(&curve),
        &MontgomeryCurve::j_invariant(&curve),
    ));
}

#[test]
fn clone_and_equality_preserve_all_coefficients() {
    let curve = MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid curve");

    let clone = curve.clone();

    assert_eq!(clone, curve);
    assert!(F7::eq(clone.a(), curve.a()));
    assert!(F7::eq(clone.b(), curve.b()));
}

#[test]
fn equation_string_mentions_every_montgomery_term() {
    let curve = MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid curve");

    let equation = curve.to_equation_string();

    assert!(equation.contains("y^2"));
    assert!(equation.contains("x^3"));
    assert!(equation.contains("x^2"));
    assert!(equation.ends_with("+ x"));
}
