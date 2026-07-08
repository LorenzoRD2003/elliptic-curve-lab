use std::collections::HashSet;

use num_bigint::BigUint;
use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    traits::{AffineCurveModel, CurveModel},
};
use crate::fields::traits::Field;
use crate::isogenies::{
    graphs::{
        IsogenyGraph, IsogenyGraphEdge, IsogenyGraphEdgeId, IsogenyGraphNode, IsogenyGraphNodeId,
        edge::EdgeTargetWitness,
        endomorphisms::{CraterShape, HorizontalEdgeStatus},
    },
    kernel::IsogenyKernel,
};
use crate::proptest_support::isogenies::arb_volcanic_floor_search_case;

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn f41_special_j_zero_curve() -> Curve41 {
    Curve41::new(F41::zero(), F41::one()).expect("valid j = 0 curve")
}

fn sample_two_torsion_kernel(curve: &Curve41) -> IsogenyKernel<Curve41> {
    let point = curve
        .point(F41::from_i64(40), F41::zero())
        .expect("sample point should lie on the curve");

    IsogenyKernel::new(curve, HashSet::from([curve.identity(), point]))
        .expect("two-torsion subgroup should be valid")
}

fn push_edge(
    edges: &mut Vec<IsogenyGraphEdge<Curve41>>,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    kernel: &IsogenyKernel<Curve41>,
) {
    edges.push(IsogenyGraphEdge::new(
        IsogenyGraphEdgeId(edges.len()),
        source,
        target,
        kernel.clone(),
        EdgeTargetWitness::Identity,
    ));
}

fn graph_from_edges(
    curve: Curve41,
    node_count: usize,
    edge_pairs: &[(usize, usize)],
    fully_expanded_nodes: Vec<bool>,
) -> IsogenyGraph<Curve41> {
    let kernel = sample_two_torsion_kernel(&curve);
    let nodes = (0..node_count)
        .map(|index| IsogenyGraphNode::new(IsogenyGraphNodeId(index), curve.clone()))
        .collect::<Vec<_>>();
    let mut edges = Vec::new();
    for &(source, target) in edge_pairs {
        push_edge(
            &mut edges,
            IsogenyGraphNodeId(source),
            IsogenyGraphNodeId(target),
            &kernel,
        );
    }

    IsogenyGraph {
        nodes,
        edges,
        fully_expanded_nodes,
    }
}

fn two_vertex_crater_graph() -> IsogenyGraph<Curve41> {
    graph_from_edges(
        f41_curve(),
        4,
        &[
            (0, 1),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 0),
            (1, 3),
            (2, 0),
            (3, 1),
        ],
        vec![true; 4],
    )
}

fn three_cycle_crater_graph() -> IsogenyGraph<Curve41> {
    graph_from_edges(
        f41_curve(),
        6,
        &[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 0),
            (1, 2),
            (1, 4),
            (2, 0),
            (2, 1),
            (2, 5),
            (3, 0),
            (4, 1),
            (5, 2),
        ],
        vec![true; 6],
    )
}

#[test]
fn crater_report_detects_two_vertex_degenerate_crater() {
    let report = two_vertex_crater_graph()
        .volcano_crater_report(&BigUint::from(2u8))
        .expect("two-vertex structural crater should report");

    assert_eq!(
        report.nodes(),
        &[IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]
    );
    assert_eq!(
        report.shape(),
        CraterShape::TwoVertex {
            directed_edge_count: 4
        }
    );
    assert_eq!(report.crater_length(), Some(2));
    assert_eq!(report.horizontal_cycle_count(), 1);
    assert_eq!(
        report.horizontal_edge_count_by_status(HorizontalEdgeStatus::CertifiedByAltitude),
        4
    );
    assert!(report.horizontal_edges().iter().all(|edge| {
        edge.status() == HorizontalEdgeStatus::CertifiedByAltitude
            && report.nodes().contains(&edge.source())
            && report.nodes().contains(&edge.target())
    }));
}

#[test]
fn crater_report_detects_simple_horizontal_cycle() {
    let report = three_cycle_crater_graph()
        .volcano_crater_report(&BigUint::from(2u8))
        .expect("three-cycle structural crater should report");

    assert_eq!(
        report.nodes(),
        &[
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2)
        ]
    );
    assert_eq!(report.shape(), CraterShape::Cycle { length: 3 });
    assert_eq!(report.crater_length(), Some(3));
    assert_eq!(
        report.horizontal_edge_count_by_status(HorizontalEdgeStatus::CertifiedByAltitude),
        6
    );
}

#[test]
fn crater_report_marks_partial_self_loop_as_not_certifiable() {
    let graph = graph_from_edges(f41_curve(), 1, &[(0, 0)], vec![false]);
    let report = graph
        .volcano_crater_report(&BigUint::from(2u8))
        .expect("valid prime should produce a crater report");

    assert_eq!(report.shape(), CraterShape::EmptyCertifiedCrater);
    assert!(report.nodes().is_empty());
    assert_eq!(report.horizontal_edges().len(), 1);
    assert_eq!(
        report.horizontal_edges()[0].status(),
        HorizontalEdgeStatus::NotCertifiableBecausePartialGraph
    );
}

#[test]
fn crater_report_marks_special_surface_self_loop_as_weakly_suspected() {
    let graph = graph_from_edges(f41_special_j_zero_curve(), 1, &[(0, 0)], vec![true]);
    let report = graph
        .volcano_crater_report(&BigUint::from(2u8))
        .expect("valid prime should produce a crater report");

    assert_eq!(report.shape(), CraterShape::EmptyCertifiedCrater);
    assert_eq!(report.horizontal_edges().len(), 1);
    assert_eq!(
        report.horizontal_edges()[0].status(),
        HorizontalEdgeStatus::SuspectedByWeakSurfaceEvidence
    );
}

proptest! {
    #[test]
    fn generated_volcanoes_have_singleton_certified_crater(
        case in arb_volcanic_floor_search_case(),
    ) {
        let report = case
            .graph()
            .volcano_crater_report(case.prime())
            .expect("complete generated volcano should produce a crater report");

        prop_assert_eq!(report.nodes(), &[IsogenyGraphNodeId(0)]);
        prop_assert_eq!(
            report.shape(),
            CraterShape::Singleton { self_loop_count: 0 }
        );
        prop_assert_eq!(report.crater_length(), Some(1));
        prop_assert_eq!(report.horizontal_cycle_count(), 1);
        prop_assert!(report.horizontal_edges().is_empty());
    }
}
