use super::shared::{F5, F13, f5_curve};
use crate::elliptic_curves::{
    AffinePoint, CurveError, TwistedEdwardsCurve,
    models::traits::{AffineCurveModel, HasProjectiveModel},
    models::twisted_edwards::projective::ExtendedTwistedEdwardsPoint,
};
use crate::fields::FieldError;
use crate::fields::traits::*;

#[test]
fn extended_identity_uses_the_expected_coordinates() {
    let identity = ExtendedTwistedEdwardsPoint::<F13>::identity();

    assert!(F13::is_zero(identity.x()));
    assert!(F13::eq(identity.y(), &F13::one()));
    assert!(F13::eq(identity.z(), &F13::one()));
    assert!(F13::is_zero(identity.t()));
    assert!(identity.is_identity());
}

#[test]
fn projective_trait_identity_matches_the_finite_twisted_edwards_neutral_element() {
    let curve = f5_curve();

    assert_eq!(
        curve.projective_identity(),
        ExtendedTwistedEdwardsPoint::<F5>::identity()
    );
}

#[test]
fn extended_from_affine_embeds_into_the_normalized_chart() {
    let point = AffinePoint::<F5>::new(F5::from_i64(2), F5::from_i64(2));
    let projective =
        ExtendedTwistedEdwardsPoint::from_affine(&point).expect("finite point should embed");

    assert!(F5::eq(projective.x(), &F5::from_i64(2)));
    assert!(F5::eq(projective.y(), &F5::from_i64(2)));
    assert!(F5::eq(projective.z(), &F5::one()));
    assert!(F5::eq(projective.t(), &F5::from_i64(4)));
}

#[test]
fn extended_from_affine_rejects_affine_infinity() {
    assert_eq!(
        ExtendedTwistedEdwardsPoint::<F13>::from_affine(&AffinePoint::Infinity),
        Err(CurveError::PointNotOnCurve)
    );
}

#[test]
fn extended_membership_accepts_affine_embeddings() {
    let curve = f5_curve();
    let affine = AffinePoint::<F5>::new(F5::from_i64(2), F5::from_i64(2));
    let extended = ExtendedTwistedEdwardsPoint::from_affine(&affine)
        .expect("finite affine point should embed");

    assert!(curve.contains_extended_point(&extended));
}

#[test]
fn extended_membership_rejects_zero_tuple() {
    let curve = f5_curve();
    let zero =
        ExtendedTwistedEdwardsPoint::<F5>::new(F5::zero(), F5::zero(), F5::zero(), F5::zero());

    assert!(!curve.contains_extended_point(&zero));
}

#[test]
fn extended_membership_rejects_structural_relation_failure() {
    let curve = f5_curve();
    let point = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(2),
        F5::from_i64(2),
        F5::one(),
        F5::zero(),
    );

    assert!(!curve.contains_extended_point(&point));
}

#[test]
fn extended_membership_rejects_curve_equation_failure() {
    let curve = f5_curve();
    let point = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(2),
        F5::from_i64(2),
        F5::one(),
        F5::from_i64(3),
    );

    assert!(!curve.contains_extended_point(&point));
}

#[test]
fn extended_curve_level_affine_roundtrip_recovers_the_original_point() {
    let curve = f5_curve();
    let affine = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");

    let extended = curve
        .extended_point_from_affine(&affine)
        .expect("affine point should lift");
    let recovered = curve
        .extended_point_to_affine(&extended)
        .expect("embedded point should recover affinely");

    assert_eq!(recovered, affine);
}

#[test]
fn projective_trait_lifts_and_recovers_affine_points() {
    let curve = f5_curve();
    let affine = curve
        .point(F5::from_i64(2), F5::from_i64(2))
        .expect("sample point should be on the curve");

    let projective = curve
        .to_projective(&affine)
        .expect("trait lift should succeed");
    let recovered = curve
        .to_affine_projective(&projective)
        .expect("trait affine recovery should succeed");

    assert_eq!(recovered, affine);
}

#[test]
fn projective_trait_reports_membership_honestly() {
    let curve = f5_curve();
    let on_curve = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(2),
        F5::from_i64(2),
        F5::one(),
        F5::from_i64(4),
    );
    let off_curve = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(2),
        F5::from_i64(2),
        F5::one(),
        F5::zero(),
    );

    assert!(curve.is_projective_point_on_curve(&on_curve));
    assert!(!curve.is_projective_point_on_curve(&off_curve));
}

#[test]
fn extended_affine_recovery_accepts_scaled_representatives() {
    let curve = f5_curve();
    let extended = ExtendedTwistedEdwardsPoint::<F5>::new(
        F5::from_i64(4),
        F5::from_i64(4),
        F5::from_i64(2),
        F5::from_i64(3),
    );

    let recovered = curve
        .extended_point_to_affine(&extended)
        .expect("scaled representative should recover affinely");

    assert_eq!(
        recovered,
        AffinePoint::new(F5::from_i64(2), F5::from_i64(2))
    );
}

#[test]
fn extended_affine_recovery_fails_honestly_when_z_is_zero() {
    let curve = TwistedEdwardsCurve::<F13>::new(F13::one(), F13::from_i64(3))
        .expect("sample twisted-Edwards curve should be non-singular");
    let point = ExtendedTwistedEdwardsPoint::<F13>::new(
        F13::zero(),
        F13::from_i64(4),
        F13::zero(),
        F13::one(),
    );

    assert!(curve.contains_extended_point(&point));
    assert_eq!(
        curve.extended_point_to_affine(&point),
        Err(CurveError::Field(FieldError::DivisionByZero))
    );
}

#[test]
fn extended_scaled_representatives_compare_equal() {
    let left = ExtendedTwistedEdwardsPoint::<F13>::new(
        F13::from_i64(2),
        F13::from_i64(5),
        F13::one(),
        F13::from_i64(10),
    );
    let right = ExtendedTwistedEdwardsPoint::<F13>::new(
        F13::from_i64(6),
        F13::from_i64(2),
        F13::from_i64(3),
        F13::from_i64(4),
    );

    assert_eq!(left, right);
    assert!(!left.has_same_representative_as(&right));
}

#[test]
fn extended_zero_tuple_only_equals_itself() {
    let zero =
        ExtendedTwistedEdwardsPoint::<F13>::new(F13::zero(), F13::zero(), F13::zero(), F13::zero());
    let identity = ExtendedTwistedEdwardsPoint::<F13>::identity();

    assert_eq!(zero, zero.clone());
    assert_ne!(zero, identity);
}

#[test]
fn extended_negation_matches_the_embedded_affine_involution() {
    let point = AffinePoint::<F13>::new(F13::from_i64(2), F13::from_i64(5));
    let negated = ExtendedTwistedEdwardsPoint::from_affine(&point)
        .expect("finite affine point should embed")
        .neg();

    let expected = ExtendedTwistedEdwardsPoint::<F13>::new(
        F13::from_i64(-2),
        F13::from_i64(5),
        F13::one(),
        F13::from_i64(-10),
    );

    assert_eq!(negated, expected);
}
