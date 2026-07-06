use num_bigint::BigUint;
use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyEdgeEndomorphismTentativeRelation, floor_search::VolcanoAltimeterEvidence,
    },
};
use crate::proptest_support::isogenies::arb_volcanic_floor_search_case;

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

proptest! {
    #[test]
    fn altimeter_distances_match_generated_volcano_levels(
        case in arb_volcanic_floor_search_case(),
    ) {
        let altimeter = VolcanoAltimeterEvidence::from_graph(case.graph(), case.prime());

        for node in case.graph().nodes() {
            let level = case
                .node_level(node.id())
                .expect("generated node should have a structural level");
            prop_assert_eq!(
                altimeter.distance_to_floor(node.id()),
                Some(case.depth() - level)
            );
        }
    }

    #[test]
    fn altimeter_classifies_generated_volcano_edges_from_delta_comparison(
        case in arb_volcanic_floor_search_case(),
    ) {
        let altimeter = VolcanoAltimeterEvidence::from_graph(case.graph(), case.prime());

        for edge in case.graph().edges() {
            let source_level = case
                .node_level(edge.source())
                .expect("generated source should have a structural level");
            let target_level = case
                .node_level(edge.target())
                .expect("generated target should have a structural level");

            prop_assert_eq!(
                altimeter.relation_for(edge.source(), edge.target()),
                IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(
                    case.depth() - source_level,
                    case.depth() - target_level,
                )
            );
        }
    }
}

#[test]
fn altimeter_leaves_partial_boundary_nodes_unclassified() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");
    let altimeter = VolcanoAltimeterEvidence::from_graph(&graph, &BigUint::from(2u8));

    assert_eq!(altimeter.distance_to_floor(IsogenyGraphNodeId(0)), None);
    assert_eq!(
        altimeter.relation_for(IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
        None
    );
}
