use crate::elliptic_curves::{
    CurveError, GeneralWeierstrassCurve,
    traits::{AffineCurveModel, CurveModel, HasJInvariant},
};
use crate::fields::traits::Field;

use super::shared::{F2, F5, F7};

#[test]
fn constructor_rejects_singular_coefficients() {
    assert!(matches!(
        GeneralWeierstrassCurve::<F5>::new(
            F5::zero(),
            F5::zero(),
            F5::zero(),
            F5::zero(),
            F5::zero(),
        ),
        Err(CurveError::SingularCurve),
    ));
}

#[test]
fn constructor_allows_characteristic_two_when_discriminant_is_nonzero() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular general Weierstrass curve in characteristic two");

    assert!(F2::eq(curve.a1(), &F2::one()));
    assert!(F2::eq(&curve.discriminant(), &F2::one()));
}

#[test]
fn short_embedding_example_matches_the_expected_classical_invariants() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    assert!(F7::eq(curve.a1(), &F7::zero()));
    assert!(F7::eq(curve.a2(), &F7::zero()));
    assert!(F7::eq(curve.a3(), &F7::zero()));
    assert!(F7::eq(curve.a4(), &F7::from_i64(2)));
    assert!(F7::eq(curve.a6(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.b2(), &F7::zero()));
    assert!(F7::eq(&curve.b4(), &F7::from_i64(4)));
    assert!(F7::eq(&curve.b6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.b8(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.c4(), &F7::from_i64(2)));
    assert!(F7::eq(&curve.c6(), &F7::from_i64(5)));
    assert!(F7::eq(&curve.discriminant(), &F7::from_i64(3)));
    assert!(F7::eq(&curve.j_invariant(), &F7::from_i64(5)));
}

#[test]
fn weierstrass_invariants_satisfy_the_classical_relation() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    let left = F7::sub(&F7::cube(&curve.c4()), &F7::square(&curve.c6()));
    let right = F7::mul(&F7::from_i64(1728), &curve.discriminant());

    assert!(F7::eq(&left, &right));
}

#[test]
fn display_and_debug_surface_the_general_equation() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let display = curve.to_string();

    assert!(display.starts_with("y^2 + ("));
    assert!(display.contains(")xy + ("));
    assert!(display.contains(")y = x^3 + ("));
    assert!(display.contains(")x^2 + ("));
    assert!(display.contains(")x + ("));
    assert!(format!("{curve:?}").contains("GeneralWeierstrassCurve"));
}

#[test]
fn curve_model_identity_and_membership_helpers_work_for_general_weierstrass() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let finite_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero());
    let infinity = crate::elliptic_curves::AffinePoint::<F5>::infinity();
    let off_curve_point = crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::one());

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
fn affine_curve_model_point_accepts_valid_general_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    let point = curve
        .point(F5::zero(), F5::zero())
        .expect("the affine point should lie on the curve");

    assert_eq!(
        point,
        crate::elliptic_curves::AffinePoint::<F5>::new(F5::zero(), F5::zero())
    );
    assert!(curve.contains(&point));
}

#[test]
fn affine_curve_model_point_rejects_invalid_general_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");

    assert_eq!(
        curve.point(F5::zero(), F5::one()),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn has_j_invariant_trait_matches_the_inherent_general_weierstrass_invariant() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::zero(),
        F7::zero(),
        F7::zero(),
        F7::from_i64(2),
        F7::from_i64(3),
    )
    .expect("non-singular curve");

    assert!(F7::eq(
        &HasJInvariant::j_invariant(&curve),
        &GeneralWeierstrassCurve::j_invariant(&curve)
    ));
}

#[test]
fn clone_and_equality_preserve_all_coefficients() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let clone = curve.clone();

    assert_eq!(clone, curve);
    assert!(F7::eq(clone.a1(), curve.a1()));
    assert!(F7::eq(clone.a2(), curve.a2()));
    assert!(F7::eq(clone.a3(), curve.a3()));
    assert!(F7::eq(clone.a4(), curve.a4()));
    assert!(F7::eq(clone.a6(), curve.a6()));
}

#[test]
fn equation_string_mentions_every_general_weierstrass_term() {
    let curve = GeneralWeierstrassCurve::<F7>::new(
        F7::one(),
        F7::from_i64(2),
        F7::from_i64(3),
        F7::from_i64(4),
        F7::from_i64(5),
    )
    .expect("non-singular curve");

    let equation = curve.to_equation_string();

    assert!(equation.contains("xy"));
    assert!(equation.contains("x^2"));
    assert!(equation.contains("x^3"));
    assert!(equation.contains("y = x^3"));
}

#[test]
fn affine_curve_model_point_accepts_valid_general_points_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    let point = curve
        .point(F2::one(), F2::zero())
        .expect("the affine point should lie on the characteristic-two curve");

    assert!(curve.contains(&point));
}
