use std::collections::HashSet;

use super::{IsogenyGraph, IsogenyGraphBuilder};
use crate::elliptic_curves::{
    AffineCurveModel, CurveIsomorphism, CurveModel, ShortWeierstrassCurve,
    ShortWeierstrassIsomorphism,
};
use crate::fields::{Field, Fp};
use crate::isogenies::graphs::{
    EdgeTargetWitness, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphError, IsogenyGraphNode,
    IsogenyGraphNodeId, ReverseEdgeStatus,
};
use crate::isogenies::{Isogeny, IsogenyKernel, VeluIsogeny};

type F7 = Fp<7>;
type F5 = Fp<5>;
type F41 = Fp<41>;
type Curve = ShortWeierstrassCurve<F41>;
type Curve5 = ShortWeierstrassCurve<F5>;

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

fn f41_curve() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn f5_split_two_torsion_curve() -> Curve5 {
    Curve5::new(F5::from_i64(-1), F5::zero()).expect("valid curve")
}

fn scaled_f41_curve() -> Curve {
    f41_curve()
        .scaled_by(F41::from_i64(3))
        .expect("non-zero scale should define a valid scaled curve")
}

fn two_torsion_kernel(curve: &Curve) -> IsogenyKernel<Curve> {
    let point = curve
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample point should lie on the curve");

    IsogenyKernel::new(curve, HashSet::from([curve.identity(), point]))
        .expect("two-torsion subgroup should be valid")
}

fn sample_graph() -> IsogenyGraph<Curve> {
    let source_curve = f41_curve();
    let target_curve = scaled_f41_curve();
    let kernel = two_torsion_kernel(&source_curve);
    let witness = ShortWeierstrassIsomorphism::new(source_curve.clone(), F41::from_i64(3))
        .expect("non-zero scale should define an isomorphism");

    IsogenyGraph {
        nodes: vec![
            IsogenyGraphNode::new(IsogenyGraphNodeId(0), source_curve),
            IsogenyGraphNode::new(IsogenyGraphNodeId(1), target_curve),
        ],
        edges: vec![
            IsogenyGraphEdge::new(
                IsogenyGraphEdgeId(0),
                IsogenyGraphNodeId(0),
                IsogenyGraphNodeId(1),
                kernel.clone(),
                EdgeTargetWitness::Explicit(witness),
            ),
            IsogenyGraphEdge::new(
                IsogenyGraphEdgeId(1),
                IsogenyGraphNodeId(1),
                IsogenyGraphNodeId(0),
                kernel,
                EdgeTargetWitness::Identity,
            ),
        ],
    }
}

fn singleton_graph() -> IsogenyGraph<Curve> {
    IsogenyGraph {
        nodes: vec![IsogenyGraphNode::new(IsogenyGraphNodeId(0), f41_curve())],
        edges: vec![],
    }
}

fn graph_with_scaled_target_only() -> IsogenyGraph<Curve> {
    IsogenyGraph {
        nodes: vec![IsogenyGraphNode::new(
            IsogenyGraphNodeId(0),
            scaled_f41_curve(),
        )],
        edges: vec![],
    }
}

fn depth_one_f41_graph() -> IsogenyGraph<Curve> {
    IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one F41 graph should build")
}

fn graph_with_actual_dual_reverse_edge() -> IsogenyGraph<Curve> {
    let mut graph = depth_one_f41_graph();
    let forward_generator = f41_curve()
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample two-torsion point should lie on the curve");
    let phi = VeluIsogeny::from_generator(f41_curve(), forward_generator)
        .expect("sample degree-two Vélu isogeny should build");
    let dual = phi
        .find_dual_exhaustively()
        .expect("small degree-two example should have an exhaustively found dual");

    graph.edges.push(IsogenyGraphEdge::new(
        IsogenyGraphEdgeId(1),
        IsogenyGraphNodeId(1),
        IsogenyGraphNodeId(0),
        dual.velu_part().kernel().clone(),
        EdgeTargetWitness::Explicit(dual.codomain_to_original().clone()),
    ));

    graph
}

fn graph_with_present_but_nondual_reverse_edge() -> IsogenyGraph<Curve> {
    let mut graph = depth_one_f41_graph();
    let forward_generator = f41_curve()
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample two-torsion point should lie on the curve");
    let phi = VeluIsogeny::from_generator(f41_curve(), forward_generator)
        .expect("sample degree-two Vélu isogeny should build");
    let dual = phi
        .find_dual_exhaustively()
        .expect("small degree-two example should have an exhaustively found dual");
    let twisted_scale = F41::mul(
        &F41::from_i64(-1),
        dual.codomain_to_original().scaling_factor(),
    );
    let nondual_witness = ShortWeierstrassIsomorphism::new(
        dual.codomain_to_original().domain().clone(),
        twisted_scale,
    )
    .expect("multiplying the dual witness by a nontrivial automorphism should stay valid");

    graph.edges.push(IsogenyGraphEdge::new(
        IsogenyGraphEdgeId(1),
        IsogenyGraphNodeId(1),
        IsogenyGraphNodeId(0),
        dual.velu_part().kernel().clone(),
        EdgeTargetWitness::Explicit(nondual_witness),
    ));

    graph
}

#[test]
fn find_isomorphic_node_returns_the_existing_node_and_witness() {
    let candidate = f41_curve();
    let representative = scaled_f41_curve();
    let nodes = vec![IsogenyGraphNode::new(
        IsogenyGraphNodeId(0),
        representative.clone(),
    )];

    let (id, witness) = IsogenyGraph::<Curve>::find_isomorphic_node(&nodes, &candidate)
        .expect("scaled curve should be found");

    assert_eq!(id, IsogenyGraphNodeId(0));
    assert_eq!(witness.domain(), &candidate);
    assert_eq!(witness.codomain(), &representative);
}

#[test]
fn find_isomorphic_node_rejects_same_j_without_base_field_isomorphism() {
    let curve = f7_curve();
    let same_j_not_base_isomorphic =
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(4), F7::from_i64(4)).expect("valid curve");
    let nodes = vec![IsogenyGraphNode::new(
        IsogenyGraphNodeId(0),
        same_j_not_base_isomorphic,
    )];

    assert!(
        IsogenyGraph::<ShortWeierstrassCurve<F7>>::find_isomorphic_node(&nodes, &curve).is_none()
    );
}

#[test]
fn resolve_target_node_reuses_existing_isomorphic_representative() {
    let mut graph = graph_with_scaled_target_only();
    let raw_codomain = f41_curve();

    let (target_id, witness) = graph.resolve_target_node(raw_codomain.clone());

    assert_eq!(target_id, IsogenyGraphNodeId(0));
    assert_eq!(graph.node_count(), 1);
    match witness {
        EdgeTargetWitness::Identity => panic!("expected explicit witness"),
        EdgeTargetWitness::Explicit(witness) => {
            assert_eq!(witness.domain(), &raw_codomain);
            assert_eq!(witness.codomain(), graph.nodes()[0].representative());
        }
    }
}

#[test]
fn graph_exposes_nodes_edges_and_dense_node_lookup() {
    let graph = sample_graph();

    assert_eq!(graph.nodes().len(), 2);
    assert_eq!(graph.node_count(), 2);
    assert_eq!(graph.edges().len(), 2);
    assert_eq!(graph.edge_count(), 2);
    assert_eq!(
        graph.node(IsogenyGraphNodeId(0)).map(|node| node.id()),
        Some(IsogenyGraphNodeId(0))
    );
    assert_eq!(
        graph.node(IsogenyGraphNodeId(1)).map(|node| node.id()),
        Some(IsogenyGraphNodeId(1))
    );
    assert!(graph.node(IsogenyGraphNodeId(2)).is_none());
}

#[test]
fn graph_stores_nodes_and_edges() {
    let graph = sample_graph();

    assert_eq!(graph.nodes().len(), 2);
    assert_eq!(graph.node_count(), 2);
    assert_eq!(graph.edges().len(), 2);
    assert_eq!(graph.edge_count(), 2);
}

#[test]
fn graph_filters_outgoing_and_incoming_edges() {
    let graph = sample_graph();

    let outgoing_from_zero = graph
        .outgoing_edges(IsogenyGraphNodeId(0))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();
    let incoming_to_zero = graph
        .incoming_edges(IsogenyGraphNodeId(0))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();

    assert_eq!(outgoing_from_zero, vec![IsogenyGraphEdgeId(0)]);
    assert_eq!(incoming_to_zero, vec![IsogenyGraphEdgeId(1)]);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(0)), 1);
    assert_eq!(graph.in_degree(IsogenyGraphNodeId(0)), 1);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(1)), 1);
    assert_eq!(graph.in_degree(IsogenyGraphNodeId(1)), 1);
}

#[test]
fn outgoing_edges_are_filtered_by_source() {
    let graph = sample_graph();

    let outgoing_from_zero = graph
        .outgoing_edges(IsogenyGraphNodeId(0))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();
    let outgoing_from_one = graph
        .outgoing_edges(IsogenyGraphNodeId(1))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();

    assert_eq!(outgoing_from_zero, vec![IsogenyGraphEdgeId(0)]);
    assert_eq!(outgoing_from_one, vec![IsogenyGraphEdgeId(1)]);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(0)), 1);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(1)), 1);
}

#[test]
fn incoming_edges_are_filtered_by_target() {
    let graph = sample_graph();

    let incoming_to_zero = graph
        .incoming_edges(IsogenyGraphNodeId(0))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();
    let incoming_to_one = graph
        .incoming_edges(IsogenyGraphNodeId(1))
        .map(|edge| edge.id())
        .collect::<Vec<_>>();

    assert_eq!(incoming_to_zero, vec![IsogenyGraphEdgeId(1)]);
    assert_eq!(incoming_to_one, vec![IsogenyGraphEdgeId(0)]);
    assert_eq!(graph.in_degree(IsogenyGraphNodeId(0)), 1);
    assert_eq!(graph.in_degree(IsogenyGraphNodeId(1)), 1);
}

#[test]
fn outgoing_velu_edges_from_node_rejects_missing_source_nodes() {
    let mut graph = singleton_graph();

    let error = graph
        .outgoing_velu_edges_from_node(IsogenyGraphNodeId(9), 2)
        .expect_err("missing source node should fail");

    assert_eq!(
        error,
        IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(9))
    );
}

#[test]
fn outgoing_velu_edges_from_node_rejects_non_prime_degrees() {
    let mut graph = singleton_graph();

    let error = graph
        .outgoing_velu_edges_from_node(IsogenyGraphNodeId(0), 4)
        .expect_err("composite degree should fail");

    assert_eq!(
        error,
        IsogenyGraphError::DegreeMustBePrimeForThisBuilder { degree: 4 }
    );
}

#[test]
fn outgoing_velu_edges_from_node_reports_missing_rational_kernels() {
    let mut graph = singleton_graph();

    let error = graph
        .outgoing_velu_edges_from_node(IsogenyGraphNodeId(0), 5)
        .expect_err("missing rational five-torsion should fail");

    assert_eq!(
        error,
        IsogenyGraphError::NonRationalKernelForRequestedDegree { degree: 5 }
    );
}

#[test]
fn outgoing_velu_edges_from_node_builds_degree_two_edges_on_f41() {
    let mut graph = singleton_graph();

    let edge_ids = graph
        .outgoing_velu_edges_from_node(IsogenyGraphNodeId(0), 2)
        .expect("degree-two rational kernel should produce one outgoing edge");

    assert_eq!(edge_ids, vec![IsogenyGraphEdgeId(0)]);
    assert_eq!(graph.node_count(), 2);
    assert_eq!(graph.edge_count(), 1);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(0)), 1);

    let edge = &graph.edges()[0];
    let expected_generator = f41_curve()
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample two-torsion point should lie on the curve");
    let expected_velu = VeluIsogeny::from_generator(f41_curve(), expected_generator)
        .expect("sample degree-two Vélu isogeny should build");

    assert_eq!(edge.source(), IsogenyGraphNodeId(0));
    assert_eq!(edge.degree(), 2);
    assert_eq!(edge.kernel(), expected_velu.kernel());
    assert_eq!(graph.nodes()[1].representative(), expected_velu.codomain());
    assert!(matches!(edge.target_witness(), EdgeTargetWitness::Identity));
}

#[test]
fn builder_build_rejects_non_prime_degree() {
    let error = IsogenyGraphBuilder::new(f41_curve(), 4)
        .build()
        .expect_err("composite degree should fail");

    assert_eq!(
        error,
        IsogenyGraphError::DegreeMustBePrimeForThisBuilder { degree: 4 }
    );
}

#[test]
fn builder_build_depth_zero_keeps_only_the_start_node() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero build should succeed");

    assert_eq!(graph.node_count(), 1);
    assert_eq!(graph.edge_count(), 0);
    assert_eq!(graph.nodes()[0].representative(), &f41_curve());
}

#[test]
fn depth_zero_graph_has_only_start_node() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero build should succeed");

    assert_eq!(graph.node_count(), 1);
    assert_eq!(graph.edge_count(), 0);
    assert_eq!(graph.nodes()[0].representative(), &f41_curve());
}

#[test]
fn builder_build_treats_missing_rational_kernels_as_leaf_behavior() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("lack of rational five-torsion should not be a global build error");

    assert_eq!(graph.node_count(), 1);
    assert_eq!(graph.edge_count(), 0);
}

#[test]
fn builder_build_depth_one_matches_the_f41_degree_two_step() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("degree-two BFS build should succeed");

    let expected_generator = f41_curve()
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample two-torsion point should lie on the curve");
    let expected_velu = VeluIsogeny::from_generator(f41_curve(), expected_generator)
        .expect("sample degree-two Vélu isogeny should build");

    assert_eq!(graph.node_count(), 2);
    assert_eq!(graph.edge_count(), 1);
    assert_eq!(graph.out_degree(IsogenyGraphNodeId(0)), 1);
    assert_eq!(graph.nodes()[1].representative(), expected_velu.codomain());
    assert_eq!(graph.edges()[0].kernel(), expected_velu.kernel());
}

#[test]
fn depth_one_graph_has_start_and_velu_codomains() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("degree-two BFS build should succeed");

    let expected_generator = f41_curve()
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("sample two-torsion point should lie on the curve");
    let expected_velu = VeluIsogeny::from_generator(f41_curve(), expected_generator)
        .expect("sample degree-two Vélu isogeny should build");

    assert_eq!(graph.node_count(), 2);
    assert_eq!(graph.edge_count(), 1);
    assert_eq!(graph.nodes()[0].representative(), &f41_curve());
    assert_eq!(graph.nodes()[1].representative(), expected_velu.codomain());
}

#[test]
fn builder_does_not_duplicate_isomorphic_codomain_nodes() {
    let deduplicated = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .build()
        .expect("depth-two BFS build should succeed");
    let without_base_field_dedup = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(false)
        .build()
        .expect("depth-two BFS build without base-field dedup should succeed");

    assert!(deduplicated.node_count() < without_base_field_dedup.node_count());
}

#[test]
fn builder_keeps_edges_to_existing_nodes() {
    let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .build()
        .expect("depth-two BFS build should succeed");

    assert!(graph.edge_count() > graph.node_count().saturating_sub(1));
}

#[test]
fn graph_builder_preserves_parallel_edges_if_kernels_are_distinct() {
    let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .build()
        .expect("depth-two BFS build should succeed");

    let parallel_edges = graph
        .outgoing_edges(IsogenyGraphNodeId(0))
        .filter(|edge| edge.target() == IsogenyGraphNodeId(1))
        .collect::<Vec<_>>();

    assert_eq!(parallel_edges.len(), 2);
    assert_ne!(parallel_edges[0].kernel(), parallel_edges[1].kernel());
}

#[test]
fn graph_contains_reverse_edges_for_known_small_example() {
    let graph = IsogenyGraphBuilder::new(f5_split_two_torsion_curve(), 2)
        .max_depth(2)
        .build()
        .expect("depth-two BFS build should succeed");

    let has_reverse_pair = graph.edges().iter().any(|edge| {
        graph.outgoing_edges(edge.target()).any(|candidate| {
            candidate.target() == edge.source() && candidate.degree() == edge.degree()
        })
    });

    assert!(has_reverse_pair);
}

#[test]
fn verify_locally_reports_missing_reverse_edges() {
    let graph = depth_one_f41_graph();

    let report = graph
        .verify_locally()
        .expect("depth-one graph should verify locally");

    assert_eq!(report.checked_edges, 1);
    assert_eq!(report.edges_mapping_domain_to_codomain, 1);
    assert_eq!(report.edges_mapping_kernel_to_identity, 1);
    assert_eq!(report.edges_homomorphism_ok, 1);
    assert_eq!(report.reverse_edges_verified_as_dual, 0);
    assert_eq!(
        report.reverse_edge_statuses,
        vec![(IsogenyGraphEdgeId(0), ReverseEdgeStatus::Missing)]
    );
}

#[test]
fn local_verification_accepts_small_depth_one_graph() {
    let graph = depth_one_f41_graph();

    let report = graph
        .verify_locally()
        .expect("depth-one graph should verify locally");

    assert_eq!(report.checked_edges, graph.edge_count());
    assert_eq!(report.edges_mapping_domain_to_codomain, graph.edge_count());
    assert_eq!(report.edges_mapping_kernel_to_identity, graph.edge_count());
    assert_eq!(report.edges_homomorphism_ok, graph.edge_count());
}

#[test]
fn verify_locally_reports_present_but_not_verified_reverse_edges() {
    let graph = graph_with_present_but_nondual_reverse_edge();

    let report = graph
        .verify_locally()
        .expect("graph with a valid but nondual reverse edge should verify locally");

    assert_eq!(report.checked_edges, 2);
    assert_eq!(report.edges_mapping_domain_to_codomain, 2);
    assert_eq!(report.edges_mapping_kernel_to_identity, 2);
    assert_eq!(report.edges_homomorphism_ok, 2);
    assert_eq!(report.reverse_edges_verified_as_dual, 0);
    assert_eq!(
        report.reverse_edge_statuses,
        vec![
            (
                IsogenyGraphEdgeId(0),
                ReverseEdgeStatus::PresentButNotVerifiedAsDual
            ),
            (
                IsogenyGraphEdgeId(1),
                ReverseEdgeStatus::PresentButNotVerifiedAsDual
            )
        ]
    );
}

#[test]
fn verify_locally_reports_verified_reverse_dual_edges() {
    let graph = graph_with_actual_dual_reverse_edge();

    let report = graph
        .verify_locally()
        .expect("graph containing an actual dual reverse edge should verify locally");

    assert_eq!(report.checked_edges, 2);
    assert_eq!(report.edges_mapping_domain_to_codomain, 2);
    assert_eq!(report.edges_mapping_kernel_to_identity, 2);
    assert_eq!(report.edges_homomorphism_ok, 2);
    assert_eq!(report.reverse_edges_verified_as_dual, 2);
    assert_eq!(
        report.reverse_edge_statuses,
        vec![
            (IsogenyGraphEdgeId(0), ReverseEdgeStatus::VerifiedAsDual),
            (IsogenyGraphEdgeId(1), ReverseEdgeStatus::VerifiedAsDual)
        ]
    );
}
