use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::CurveIsomorphismError,
    traits::{CurveIsomorphism, CurveModel, EnumerableCurveModel},
};
use crate::fields::traits::SqrtField;
use crate::fields::traits::*;

use super::shared::{
    F7, F13, F19, F37, f7_curve, f7_point, f13_generic_curve, f13_j_1728_curve, f13_j_zero_curve,
    f19_curve, f37_curve, first_nonsquare,
};

#[test]
fn scaled_by_rejects_noninvertible_scale() {
    let curve = f7_curve();

    assert!(matches!(
        curve.scaled_by(F7::zero()),
        Err(CurveIsomorphismError::NonInvertibleScale)
    ));
}

#[test]
fn scaled_by_applies_the_expected_u4_and_u6_coefficients() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F7::eq(scaled.a(), &F7::from_i64(1)));
    assert!(F7::eq(scaled.b(), &F7::from_i64(3)));
}

#[test]
fn isomorphic_via_scale_matches_the_scaled_curve() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.isomorphic_via_scale(&scaled, &F7::from_i64(3)));
    assert!(!curve.isomorphic_via_scale(&scaled, &F7::from_i64(2)));
}

#[test]
fn isomorphic_via_scale_returns_false_for_noninvertible_scale() {
    let curve = f7_curve();

    assert!(!curve.isomorphic_via_scale(&curve, &F7::zero()));
}

#[test]
fn quadratic_twist_rejects_zero_factor() {
    let curve = f19_curve();

    assert!(matches!(
        curve.quadratic_twist(F19::zero()),
        Err(CurveIsomorphismError::NonInvertibleScale)
    ));
}

#[test]
fn quadratic_twist_preserves_the_j_invariant_over_f19() {
    let curve = f19_curve();
    let twist = curve
        .quadratic_twist(F19::from_i64(2))
        .expect("non-zero twist factor should produce a valid model");

    assert!(curve.has_same_j_invariant(&twist));
}

#[test]
fn quadratic_twist_by_a_square_is_base_field_isomorphic_over_f19() {
    let curve = f19_curve();
    let square = F19::from_i64(4);
    let twist = curve
        .quadratic_twist(square)
        .expect("square twist factor should produce a valid model");

    assert!(F19::has_square_root(&square));
    assert!(curve.is_isomorphic_to(&twist));
}

#[test]
fn quadratic_twist_by_a_nonsquare_is_not_base_field_isomorphic_in_the_sample_f19_case() {
    let curve = f19_curve();
    let nonsquare = first_nonsquare::<F19>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert!(!F19::has_square_root(&nonsquare));
    assert!(curve.has_same_j_invariant(&twist));
    assert!(!curve.is_isomorphic_to(&twist));
}

#[test]
fn quadratic_twist_point_count_relation_holds_over_f19() {
    let curve = f19_curve();
    let nonsquare = first_nonsquare::<F19>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert_eq!(curve.order() + twist.order(), 2 * 19 + 2);
}

#[test]
fn quadratic_twist_point_count_relation_holds_over_f37() {
    let curve = f37_curve();
    let nonsquare = first_nonsquare::<F37>();
    let twist = curve
        .quadratic_twist(nonsquare)
        .expect("non-zero twist factor should produce a valid model");

    assert_eq!(curve.order() + twist.order(), 2 * 37 + 2);
}

#[test]
fn find_isomorphism_to_recovers_a_base_field_scaling_witness() {
    let curve = f7_curve();
    let other = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");
    let point = f7_point(2, 1);

    let isomorphism = curve
        .find_isomorphism_to(&other)
        .expect("a base-field scaling witness should exist");
    let image = isomorphism
        .evaluate(&point)
        .expect("the witness isomorphism should transport domain points");

    assert!(isomorphism.codomain().contains(&image));
    assert_eq!(
        isomorphism
            .inverse()
            .expect("inverse should exist")
            .evaluate(&image)
            .expect("inverse should recover the original point"),
        point
    );
}

#[test]
fn find_isomorphism_to_returns_none_for_same_j_but_no_base_field_scale() {
    let curve = f7_curve();
    let same_j_not_base_isomorphic =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4)).expect("valid curve");

    assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
    assert!(
        curve
            .find_isomorphism_to(&same_j_not_base_isomorphic)
            .is_none()
    );
}

#[test]
fn is_isomorphic_to_returns_true_when_a_base_field_witness_exists() {
    let curve = f7_curve();
    let other = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.is_isomorphic_to(&other));
}

#[test]
fn is_isomorphic_to_returns_false_when_only_the_j_invariant_matches() {
    let curve = f7_curve();
    let same_j_not_base_isomorphic =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4)).expect("valid curve");

    assert!(curve.has_same_j_invariant(&same_j_not_base_isomorphic));
    assert!(!curve.is_isomorphic_to(&same_j_not_base_isomorphic));
}

#[test]
fn generic_curve_has_only_plus_minus_one_automorphisms_over_f7() {
    let curve = f7_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 2);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F7::eq(iso.scaling_factor(), &F7::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F7::eq(iso.scaling_factor(), &F7::from_i64(-1)))
    );
}

#[test]
fn special_j_1728_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_j_1728_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F13::is_zero(curve.b()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_1728_curve_has_four_automorphisms_over_f13() {
    let curve = f13_j_1728_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 4);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(5)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-5)))
    );
}

#[test]
fn special_j_zero_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_j_zero_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(F13::is_zero(curve.a()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_zero_curve_has_six_automorphisms_over_f13() {
    let curve = f13_j_zero_curve();
    let automorphisms = curve.automorphisms();

    assert_eq!(automorphisms.len(), 6);
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::one()))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-1)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(4)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-4)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(3)))
    );
    assert!(
        automorphisms
            .iter()
            .any(|iso| F13::eq(iso.scaling_factor(), &F13::from_i64(-3)))
    );
}

#[test]
fn generic_family_supports_exhaustive_base_field_isomorphism_search() {
    let curve = f13_generic_curve();
    let scaled = curve
        .scaled_by(F13::from_i64(2))
        .expect("non-zero scale should define a valid scaled model");

    assert!(!F13::is_zero(curve.a()));
    assert!(!F13::is_zero(curve.b()));
    assert!(curve.has_same_j_invariant(&scaled));
    assert!(curve.is_isomorphic_to(&scaled));
    assert!(curve.find_isomorphism_to(&scaled).is_some());
}

#[test]
fn special_j_1728_same_j_does_not_force_base_field_isomorphism() {
    let first = f13_j_1728_curve();
    let second =
        ShortWeierstrassCurve::<F13>::new(F13::from_i64(1), F13::zero()).expect("valid curve");

    assert!(first.has_same_j_invariant(&second));
    assert!(!first.is_isomorphic_to(&second));
    assert!(first.find_isomorphism_to(&second).is_none());
}

#[test]
fn special_j_zero_same_j_does_not_force_base_field_isomorphism() {
    let first = f13_j_zero_curve();
    let second =
        ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(1)).expect("valid curve");

    assert!(first.has_same_j_invariant(&second));
    assert!(!first.is_isomorphic_to(&second));
    assert!(first.find_isomorphism_to(&second).is_none());
}

#[test]
fn has_same_j_invariant_detects_scaled_models() {
    let curve = f7_curve();
    let scaled = curve
        .scaled_by(F7::from_i64(3))
        .expect("non-zero scale should define a valid scaled model");

    assert!(curve.has_same_j_invariant(&scaled));
}

#[test]
fn has_same_j_invariant_distinguishes_curves_with_different_j() {
    let first = f7_curve();
    let second =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::from_i64(1)).expect("valid curve");

    assert!(!first.has_same_j_invariant(&second));
}
