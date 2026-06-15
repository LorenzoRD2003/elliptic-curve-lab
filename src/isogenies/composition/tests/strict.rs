use crate::elliptic_curves::{AffinePoint, CurveError};
use crate::fields::traits::Field;
use crate::isogenies::{
    comparison::maps_equal_exhaustively,
    composition::ComposedIsogeny,
    error::{IsogenyError, IsogenyMapError, IsogenyVerificationError},
    scalar_multiplication::ScalarMultiplicationIsogeny,
    traits::{Isogeny, VerifiableIsogeny},
};

use super::shared::{
    BrokenMiddleImageIsogeny, F41, curve_a, first_isogeny, first_non_kernel_point, second_isogeny,
};

#[test]
fn new_strict_accepts_exact_matching_middle_curves() {
    let first = first_isogeny();
    let second = second_isogeny(first.codomain());

    let composed = ComposedIsogeny::new_strict(first, second).expect("composition should build");

    assert_eq!(composed.first().codomain(), composed.second().domain());
}

#[test]
fn new_strict_rejects_mismatched_middle_curves() {
    let first = first_isogeny();
    let unrelated_second = second_isogeny(&curve_a());

    assert!(matches!(
        ComposedIsogeny::new_strict(first, unrelated_second),
        Err(IsogenyError::Map(
            IsogenyMapError::CompositionDomainCodomainMismatch
        ))
    ));
}

#[test]
fn degree_multiplies_the_component_degrees() {
    let first = first_isogeny();
    let second = second_isogeny(first.codomain());
    let expected_degree = first.degree() * second.degree();

    let composed = ComposedIsogeny::new_strict(first, second).expect("composition should build");

    assert_eq!(composed.degree(), expected_degree);
}

#[test]
fn evaluate_applies_first_then_second_on_the_shared_middle_curve() {
    let first = first_isogeny();
    let second = second_isogeny(first.codomain());
    let point = first_non_kernel_point(&first);
    let expected = second
        .evaluate(
            &first
                .evaluate(&point)
                .expect("first Vélu isogeny should evaluate on sample point"),
        )
        .expect("second Vélu isogeny should evaluate on the transported point");

    let composed = ComposedIsogeny::new_strict(first, second).expect("composition should build");

    assert_eq!(
        composed
            .evaluate(&point)
            .expect("composition should evaluate"),
        expected
    );
}

#[test]
fn evaluate_rejects_points_outside_the_first_domain_before_delegating() {
    let first = first_isogeny();
    let second = second_isogeny(first.codomain());
    let invalid = AffinePoint::new(F41::from_i64(2), F41::from_i64(2));

    let composed = ComposedIsogeny::new_strict(first, second).expect("composition should build");

    assert_eq!(
        composed.evaluate(&invalid),
        Err(IsogenyError::Curve(CurveError::PointNotOnCurve))
    );
}

#[test]
fn constructor_rejects_intermediate_images_outside_the_shared_middle_curve() {
    let inner_first = first_isogeny();
    let point = first_non_kernel_point(&inner_first);
    let broken_first = BrokenMiddleImageIsogeny {
        inner: inner_first,
        broken_point: point.clone(),
    };
    let second = second_isogeny(broken_first.codomain());

    assert!(matches!(
        ComposedIsogeny::new_strict(broken_first, second),
        Err(IsogenyError::Verification(
            IsogenyVerificationError::ImagePointNotOnCodomain
        ))
    ));
}

#[test]
fn strict_composition_kernel_points_match_the_full_identity_fiber() {
    let first = first_isogeny();
    let second = second_isogeny(first.codomain());
    let composed = ComposedIsogeny::new_strict(first, second).expect("composition should build");

    assert_eq!(composed.verify_kernel_exactness(), Ok(()));
}

#[test]
fn composition_with_identity_like_scalar_one_is_neutral() {
    let first = first_isogeny();
    let identity_like = ScalarMultiplicationIsogeny::new(first.codomain().clone(), 1)
        .expect("multiplication by one should build");
    let composed = ComposedIsogeny::new_strict(first.clone(), identity_like)
        .expect("composition with scalar-one should build");

    assert_eq!(
        maps_equal_exhaustively::<_, _, _, _>(&composed, &first),
        Ok(true)
    );
}
