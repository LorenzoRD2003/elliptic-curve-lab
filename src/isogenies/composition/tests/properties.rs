use proptest::prelude::*;

use crate::elliptic_curves::traits::CurveIsomorphism;
use crate::isogenies::{composition::ComposedIsogeny, traits::Isogeny};
use crate::proptest_support::isogenies::arb_composable_velu_case;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn property_strict_composition_matches_manual_sequential_evaluation(
        case in arb_composable_velu_case(),
    ) {
        let composed = ComposedIsogeny::new_strict(
            case.first.clone(),
            case.second_strict.clone(),
        )
            .expect("composition should build");
        let expected = case.second_strict
            .evaluate(
                &case.first
                    .evaluate(&case.sample_point)
                    .expect("first isogeny should evaluate"),
            )
            .expect("second isogeny should evaluate");

        prop_assert_eq!(
            composed
                .evaluate(&case.sample_point)
                .expect("composed isogeny should evaluate"),
            expected
        );
    }

    #[test]
    fn property_bridged_composition_matches_manual_bridge_transport(
        case in arb_composable_velu_case(),
    ) {
        let composed = ComposedIsogeny::new_up_to_isomorphism(
            case.first.clone(),
            case.bridge.clone(),
            case.second_bridged.clone(),
        )
        .expect("bridged composition should build");
        let expected = case.second_bridged
            .evaluate(
                &case.bridge
                    .evaluate(
                        &case.first
                            .evaluate(&case.sample_point)
                            .expect("first isogeny should evaluate"),
                    )
                    .expect("bridge should evaluate"),
            )
            .expect("second isogeny should evaluate");

        prop_assert_eq!(
            composed
                .evaluate(&case.sample_point)
                .expect("composed isogeny should evaluate"),
            expected
        );
    }
}
