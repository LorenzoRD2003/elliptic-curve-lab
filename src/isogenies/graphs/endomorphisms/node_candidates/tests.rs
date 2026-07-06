use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphError, IsogenyGraphNodeId};

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

#[test]
fn node_can_derive_its_candidate_endomorphism_orders() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");
    let node = graph
        .node(IsogenyGraphNodeId(0))
        .expect("root node should exist");

    let candidate_set = node
        .endomorphism_ring_candidates()
        .expect("candidate endomorphism orders should derive");

    assert!(!candidate_set.is_empty());
    assert_eq!(
        candidate_set.maximal_order().conductor(),
        &num_bigint::BigUint::from(1u8)
    );
}

#[test]
fn graph_can_derive_candidates_for_one_requested_node() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let candidate_set = graph
        .node_endomorphism_candidates(IsogenyGraphNodeId(0))
        .expect("root candidate set should derive");

    assert!(!candidate_set.is_empty());
}

#[test]
fn graph_can_derive_candidates_for_every_stored_node() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let candidate_sets = graph
        .graph_endomorphism_candidates()
        .expect("all node candidate sets should derive");

    assert_eq!(candidate_sets.len(), graph.node_count());
    assert_eq!(candidate_sets[0].0, IsogenyGraphNodeId(0));
    assert!(
        candidate_sets
            .iter()
            .all(|(_, candidate_set)| !candidate_set.is_empty())
    );
}

#[test]
fn missing_node_reports_a_graph_error() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");

    assert_eq!(
        graph.node_endomorphism_candidates(IsogenyGraphNodeId(99)),
        Err(IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(99)))
    );
}
