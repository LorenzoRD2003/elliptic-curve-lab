use num_bigint::BigUint;
use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId,
    endomorphisms::{VolcanoSearchError, VolcanoStructureRole, VolcanoStructureUncertifiedReason},
};
use crate::numerics::PositivePrimeError;
use crate::proptest_support::isogenies::arb_volcanic_floor_search_case;

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

#[test]
fn structure_report_separates_partial_nodes_from_certified_levels() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");

    let report = graph
        .volcano_structure_report(&BigUint::from(2u8))
        .expect("structural report should build for a valid prime");

    assert_eq!(report.prime(), &BigUint::from(2u8));
    assert_eq!(report.certified_depth(), None);
    assert!(report.levels().is_empty());
    assert!(report.certified_nodes().is_empty());
    assert!(!report.is_fully_certified());
    assert_eq!(report.uncertified_nodes().len(), graph.node_count());
    assert_eq!(
        report
            .uncertified_node(IsogenyGraphNodeId(0))
            .expect("root should be uncertified")
            .reason(),
        &VolcanoStructureUncertifiedReason::PartialGraph
    );
}

#[test]
fn structure_report_handles_depth_zero_volcanoes_honestly() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");

    let report = graph
        .volcano_structure_report(&BigUint::from(5u8))
        .expect("structural report should build");

    assert_eq!(report.certified_depth(), Some(0));
    assert_eq!(report.surface(), report.floor());
    assert_eq!(report.crater(), report.surface());
    assert_eq!(
        report
            .surface()
            .expect("depth-zero volcano has one level")
            .role(),
        VolcanoStructureRole::SurfaceAndFloor
    );
    assert_eq!(report.certified_nodes().len(), 1);
    assert_eq!(
        report.certified_nodes()[0].role(),
        VolcanoStructureRole::SurfaceAndFloor
    );
    assert!(report.is_fully_certified());
}

#[test]
fn structure_report_rejects_non_prime_local_parameter() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("graph should build");

    let error = graph
        .volcano_structure_report(&BigUint::from(4u8))
        .expect_err("composite local parameter should fail");

    assert_eq!(
        error,
        VolcanoSearchError::InvalidLocalPrime(PositivePrimeError::Composite)
    );
}

proptest! {
    #[test]
    fn structure_report_recovers_generated_volcano_levels(
        case in arb_volcanic_floor_search_case(),
    ) {
        let report = case
            .graph()
            .volcano_structure_report(case.prime())
            .expect("complete generated volcano should report structurally");

        prop_assert_eq!(report.certified_depth(), Some(case.depth()));
        prop_assert_eq!(report.levels().len(), case.depth() + 1);
        prop_assert_eq!(report.certified_nodes().len(), case.graph().node_count());
        prop_assert!(report.uncertified_nodes().is_empty());
        prop_assert!(report.is_fully_certified());

        for node in report.certified_nodes() {
            let expected_level = case
                .node_level(node.node_id())
                .expect("generated node should have a structural level");

            prop_assert_eq!(node.level(), expected_level);
            prop_assert_eq!(node.distance_to_floor(), case.depth() - expected_level);
            prop_assert_eq!(node.floor_path().distance_to_floor(), node.distance_to_floor());

            let expected_role = if case.depth() == 0 {
                VolcanoStructureRole::SurfaceAndFloor
            } else if expected_level == 0 {
                VolcanoStructureRole::Surface
            } else if expected_level == case.depth() {
                VolcanoStructureRole::Floor
            } else {
                VolcanoStructureRole::Middle
            };
            prop_assert_eq!(node.role(), expected_role);
        }

        for level in report.levels() {
            prop_assert_eq!(
                level.nodes().iter().all(|node_id| {
                    case.node_level(*node_id) == Some(level.level())
                }),
                true
            );
        }

        let floor_nodes = report
            .floor()
            .expect("complete generated volcano has a certified floor")
            .nodes();
        prop_assert_eq!(floor_nodes, case.floor_nodes());
    }
}
