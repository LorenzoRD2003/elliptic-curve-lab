use crate::elliptic_curves::traits::CurveIsomorphism;
use crate::isogenies::{
    composition::ComposedIsogeny,
    error::{IsogenyError, IsogenyMapError},
    traits::Isogeny,
};

use super::shared::{bridged_second_isogeny, curve_a, first_isogeny, first_non_kernel_point};
use crate::elliptic_curves::short_weierstrass::isomorphisms::ShortWeierstrassIsomorphism;

#[test]
fn new_up_to_isomorphism_accepts_a_valid_bridge_and_evaluates_as_psi_alpha_phi() {
    let first = first_isogeny();
    let point = first_non_kernel_point(&first);
    let middle_image = first
        .evaluate(&point)
        .expect("first Vélu isogeny should evaluate on sample point");
    let (bridge, second) = bridged_second_isogeny(&first);
    let expected = second
        .evaluate(
            &bridge
                .evaluate(&middle_image)
                .expect("bridge should transport the middle image"),
        )
        .expect("second Vélu isogeny should evaluate on the bridged point");

    let composed = ComposedIsogeny::new_up_to_isomorphism(first, bridge, second)
        .expect("bridged composition should build");

    assert_eq!(
        composed
            .evaluate(&point)
            .expect("bridged composition should evaluate"),
        expected
    );
}

#[test]
fn new_up_to_isomorphism_rejects_a_bridge_with_the_wrong_domain() {
    let first = first_isogeny();
    let wrong_bridge = ShortWeierstrassIsomorphism::new(curve_a(), super::shared::F41::from_i64(3))
        .expect("sample wrong bridge should still be a valid isomorphism");
    let second = super::shared::second_isogeny(wrong_bridge.codomain());

    assert!(matches!(
        ComposedIsogeny::new_up_to_isomorphism(first, wrong_bridge, second),
        Err(IsogenyError::Map(
            IsogenyMapError::CompositionDomainCodomainMismatch
        ))
    ));
}
