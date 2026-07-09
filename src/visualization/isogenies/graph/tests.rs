use std::collections::HashSet;

use num_bigint::BigUint;

use super::*;
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::Field;
use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};

type F5 = crate::fields::Fp5;
type F41 = crate::fields::Fp41;
type Curve5 = ShortWeierstrassCurve<F5>;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn f5_split_two_torsion_curve() -> Curve5 {
    Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
}

#[test]
fn summary_reports_depth_zero_graph_shape() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");

    let summary = graph.summary();

    assert_eq!(summary.node_count, 1);
    assert_eq!(summary.edge_count, 0);
    assert_eq!(summary.degree, 0);
    assert_eq!(summary.connected_component_count, 1);
    assert!(!summary.has_directed_cycle);
    assert_eq!(summary.self_loops, 0);
    assert_eq!(summary.repeated_j_invariants, 0);
    assert_eq!(summary.min_out_degree, 0);
    assert_eq!(summary.max_out_degree, 0);
    assert_eq!(summary.volcano_like.root, Some(IsogenyGraphNodeId(0)));
    assert_eq!(
        summary.volcano_like.levels(),
        vec![vec![IsogenyGraphNodeId(0)]]
    );
    assert_eq!(summary.volcano_like.surface_nodes, 0);
    assert_eq!(summary.volcano_like.middle_nodes, 0);
    assert_eq!(summary.volcano_like.floor_nodes, 0);
    assert_eq!(summary.volcano_like.isolated_nodes, 1);
    assert_eq!(summary.volcano_like.unknown_nodes, 0);
}

#[test]
fn summary_reports_depth_one_f41_graph_shape() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let summary = graph.summary();

    assert_eq!(
        (summary.node_count, summary.edge_count, summary.degree),
        (2, 1, 2)
    );
    assert_eq!(summary.connected_component_count, 1);
    assert!(!summary.has_directed_cycle);
    assert_eq!(summary.self_loops, 0);
    assert_eq!(summary.repeated_j_invariants, 0);
    assert_eq!(summary.min_out_degree, 0);
    assert_eq!(summary.max_out_degree, 1);
    assert_eq!(summary.volcano_like.root, Some(IsogenyGraphNodeId(0)));
    assert_eq!(
        summary.volcano_like.levels(),
        vec![vec![IsogenyGraphNodeId(0)], vec![IsogenyGraphNodeId(1)]]
    );
    assert_eq!(summary.volcano_like.surface_nodes, 0);
    assert_eq!(summary.volcano_like.middle_nodes, 0);
    assert_eq!(summary.volcano_like.floor_nodes, 2);
    assert_eq!(summary.volcano_like.isolated_nodes, 0);
    assert_eq!(summary.volcano_like.unknown_nodes, 0);
}

#[test]
fn summary_detects_repeated_j_invariants_in_split_two_torsion_example() {
    let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(false)
        .build()
        .expect("split two-torsion graph should build");

    let summary = graph.summary();
    let unique_j_count = graph
        .nodes()
        .iter()
        .map(|node| node.j_invariant())
        .collect::<HashSet<_>>()
        .len();

    assert_eq!(
        summary.repeated_j_invariants,
        summary.node_count.saturating_sub(unique_j_count)
    );
    assert_eq!(summary.connected_component_count, 1);
    assert!(summary.has_directed_cycle);
    assert_eq!(summary.degree, 2);
    assert!(summary.volcano_like.root.is_some());
    assert!(!summary.volcano_like.is_empty());
}

#[test]
fn adjacency_list_formats_dense_node_order() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let adjacency = format_adjacency_list(&graph);

    assert!(adjacency.contains("v0 -> v1"));
    assert!(adjacency.contains("v1 ->"));
}

#[test]
fn graph_explanation_mentions_summary_nodes_edges_and_adjacency() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let explanation = explain_isogeny_graph(&graph);

    assert!(explanation.contains("ℓ-isogeny graph summary"));
    assert!(explanation.contains("degree ℓ: 2"));
    assert!(explanation.contains("nodes: 2"));
    assert!(explanation.contains("edges: 1"));
    assert!(explanation.contains("has directed cycle: no"));
    assert!(explanation.contains("volcano-like root: v0"));
    assert!(explanation.contains("volcano-like levels: 2"));
    assert!(explanation.contains("Nodes:"));
    assert!(explanation.contains("Edges:"));
    assert!(explanation.contains("Adjacency list:"));
    assert!(explanation.contains("Volcano-like levels (heuristic):"));
    assert!(explanation.contains("v0: j = "));
    assert!(explanation.contains("curve = y^2 = x^3"));
    assert!(explanation.contains("e0: v0 -> v1, degree 2, kernel size 2"));
}

#[test]
fn graph_verification_explanation_summarizes_reverse_edge_statuses() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");
    let report = graph
        .verify_locally()
        .expect("tiny graph verification should run");

    let explanation = explain_graph_verification_report(&report);

    assert!(explanation.contains("Local graph verification report"));
    assert!(explanation.contains("checked edges: 1"));
    assert!(explanation.contains("maps domain to codomain: 1/1"));
    assert!(explanation.contains("maps kernel to identity: 1/1"));
    assert!(explanation.contains("homomorphism law verified: 1/1"));
    assert!(
        explanation
            .contains("reverse-edge statuses: verified 0, present-not-verified 0, missing 1")
    );
    assert!(!explanation.contains("Reverse-edge details:"));
}

#[test]
fn volcano_layering_explanation_mentions_levels_and_roles() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");
    let layers = graph.infer_volcano_like_layers(IsogenyGraphNodeId(0));

    let explanation = explain_volcano_like_layers(&graph, &layers);

    assert!(explanation.contains("Volcano-like layering (heuristic)"));
    assert!(explanation.contains("levels: 2"));
    assert!(explanation.contains("Levels:"));
    assert!(explanation.contains("Node roles:"));
    assert!(explanation.contains("v0: Floor"));
    assert!(explanation.contains("v1: Floor"));
}

#[test]
fn graph_endomorphism_report_explanation_mentions_tentative_arithmetic_data() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");
    let report = graph
        .endomorphism_report_at(&BigUint::from(2u8))
        .expect("endomorphism report should build");

    let explanation = explain_graph_endomorphism_report(&report);

    assert!(explanation.contains("Tentative endomorphism-side report"));
    assert!(explanation.contains("prime ℓ: 2"));
    assert!(explanation.contains("Frobenius-compatible only"));
    assert!(explanation.contains("Nodes:"));
    assert!(explanation.contains("possible levels"));
    assert!(explanation.contains("Edges:"));
    assert!(explanation.contains("source levels"));
    assert!(explanation.contains("target levels"));
}

#[test]
fn summary_type_is_cloneable_and_debuggable() {
    let summary = IsogenyGraphSummary {
        node_count: 1,
        edge_count: 2,
        degree: 3,
        connected_component_count: 4,
        has_directed_cycle: true,
        self_loops: 5,
        repeated_j_invariants: 6,
        min_out_degree: 7,
        max_out_degree: 8,
        volcano_like: VolcanoHeuristicSummary {
            root: None,
            levels: Vec::new(),
            surface_nodes: 0,
            middle_nodes: 0,
            floor_nodes: 0,
            isolated_nodes: 0,
            unknown_nodes: 0,
        },
    };

    let clone = summary.clone();
    let debug = format!("{summary:?}");

    assert_eq!(clone.node_count, 1);
    assert!(debug.contains("connected_component_count"));
    assert!(debug.contains("VolcanoHeuristicSummary"));
}
