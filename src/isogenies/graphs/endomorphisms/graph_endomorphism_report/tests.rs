use num_bigint::BigUint;
use proptest::prelude::*;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyEdgeEndomorphismTentativeRelation, IsogenyGraphEndomorphismReport,
        refinement::{
            CandidateRefinementError, CandidateRefinementStrategy, ConstraintSource,
            EndomorphismCandidateRefinement, IsogenyGraphCandidateRefinementReport,
            LocalEndomorphismConstraint, RefinementConfidence,
        },
    },
};

type F41 = crate::fields::Fp41;
type F17 = crate::fields::Fp17;
type Curve41 = ShortWeierstrassCurve<F41>;
type Curve17 = ShortWeierstrassCurve<F17>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

fn f17_curve(a: i64, b: i64) -> Curve17 {
    Curve17::new(F17::from_i64(a), F17::from_i64(b)).expect("valid curve")
}

fn graph_report(
    depth: usize,
) -> (
    crate::isogenies::graphs::IsogenyGraph<Curve41>,
    IsogenyGraphEndomorphismReport,
) {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(depth)
        .build()
        .expect("graph should build from the concrete curve");
    let report = graph
        .endomorphism_report_at(&BigUint::from(2u8))
        .expect("graph endomorphism report should build");
    (graph, report)
}

fn f17_graph_report(
    a: i64,
    b: i64,
    depth: usize,
) -> (
    crate::isogenies::graphs::IsogenyGraph<Curve17>,
    IsogenyGraphEndomorphismReport,
) {
    let graph = IsogenyGraphBuilder::new(f17_curve(a, b), 2)
        .max_depth(depth)
        .build()
        .expect("F17 graph should build");
    let report = graph
        .endomorphism_report_at(&BigUint::from(2u8))
        .expect("F17 graph endomorphism report should build");
    (graph, report)
}

fn refinement_at(
    depth: usize,
    strategy: CandidateRefinementStrategy,
) -> EndomorphismCandidateRefinement {
    let (_, report) = graph_report(depth);
    report
        .refine_candidates_for_node(IsogenyGraphNodeId(0), strategy)
        .expect("candidate refinement should build")
}

fn assert_only_node_or_unambiguous_edge_constraints(refinement: &EndomorphismCandidateRefinement) {
    assert!(refinement.constraints().iter().all(|constraint| {
        matches!(constraint, LocalEndomorphismConstraint::NodeLevel { .. })
            || matches!(
                constraint,
                LocalEndomorphismConstraint::EdgeRelation {
                    relation: IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
                        | IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending
                        | IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending,
                    ..
                }
            )
    }));
}

fn assert_refinement_candidates_stay_inside_initial_set(
    refinement: &EndomorphismCandidateRefinement,
) {
    let initial_orders = refinement.initial_candidates().candidate_orders();
    assert!(!initial_orders.is_empty());
    assert!(
        refinement
            .surviving_candidates()
            .iter()
            .all(|candidate| initial_orders.contains(candidate))
    );
    assert!(
        refinement
            .eliminated_candidates()
            .iter()
            .all(|elimination| initial_orders.contains(elimination.candidate()))
    );
}

fn total_initial_candidates(report: &IsogenyGraphCandidateRefinementReport) -> usize {
    report
        .node_refinements()
        .iter()
        .map(|node| node.initial_candidates().candidate_orders().len())
        .sum()
}

fn total_surviving_candidates(report: &IsogenyGraphCandidateRefinementReport) -> usize {
    report
        .node_refinements()
        .iter()
        .map(|node| node.surviving_candidates().len())
        .sum()
}

fn assert_fixed_point_is_subset_of(
    fixed_point: &IsogenyGraphCandidateRefinementReport,
    independent: &IsogenyGraphCandidateRefinementReport,
) {
    assert_eq!(
        fixed_point.node_refinements().len(),
        independent.node_refinements().len()
    );
    for fixed_node in fixed_point.node_refinements() {
        let independent_node = independent
            .refinement_for_node(fixed_node.node_id())
            .expect("matching independent node refinement should exist");
        assert!(
            fixed_node
                .surviving_candidates()
                .iter()
                .all(|candidate| independent_node.surviving_candidates().contains(candidate))
        );
    }
}

#[test]
fn graph_report_collects_node_and_edge_endomorphism_data() {
    let (graph, report) = graph_report(1);

    assert_eq!(report.prime(), &BigUint::from(2u8));
    assert_eq!(report.nodes().len(), graph.node_count());
    assert_eq!(report.edges().len(), graph.edge_count());
    assert!(
        report
            .nodes()
            .iter()
            .all(|node| !node.candidate_set().is_empty())
    );
    assert!(
        report.nodes().iter().all(
            |node| node.local_level_candidate_count() > 0 && !node.possible_levels().is_empty()
        )
    );
}

#[test]
fn graph_report_exposes_dense_node_lookup() {
    let (_, report) = graph_report(0);

    let node_report = report
        .node_report(IsogenyGraphNodeId(0))
        .expect("root node report should exist");

    assert_eq!(node_report.node_id(), IsogenyGraphNodeId(0));
    assert_eq!(
        node_report.candidate_set(),
        report.nodes()[0].candidate_set()
    );
    assert!(report.node_report(IsogenyGraphNodeId(99)).is_none());
}

#[test]
fn graph_report_edge_endomorphism_relations_are_tentative() {
    let (_, report) = graph_report(1);

    assert!(report.edges().iter().all(|edge| {
        matches!(
            edge.relation().relation(),
            IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
                | IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending
                | IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending
                | IsogenyEdgeEndomorphismTentativeRelation::Ambiguous
                | IsogenyEdgeEndomorphismTentativeRelation::Unsupported
        )
    }));
}

#[test]
fn node_local_refinement_keeps_candidates_supported_by_node_levels() {
    let (_, report) = graph_report(0);
    let node_report = report
        .node_report(IsogenyGraphNodeId(0))
        .expect("root node report should exist");

    let refinement = refinement_at(0, CandidateRefinementStrategy::NodeLocalLevelsOnly);

    assert_eq!(refinement.node_id(), IsogenyGraphNodeId(0));
    assert_eq!(refinement.initial_candidates(), node_report.candidate_set());
    assert_eq!(
        refinement.surviving_candidates(),
        node_report.candidate_set().candidate_orders()
    );
    assert!(refinement.eliminated_candidates().is_empty());
    assert_eq!(
        refinement.confidence(),
        RefinementConfidence::ConservativeLocalEvidence
    );
    assert_eq!(refinement.constraints().len(), 1);

    let LocalEndomorphismConstraint::NodeLevel {
        ell,
        allowed_levels,
        provenance,
    } = &refinement.constraints()[0]
    else {
        panic!("node-local refinement should record one node-level constraint");
    };
    let expected_levels = node_report.possible_levels().into_iter().collect();
    assert_eq!(ell, &BigUint::from(2u8));
    assert_eq!(allowed_levels, &expected_levels);
    assert_eq!(
        provenance,
        &ConstraintSource::NodeReport {
            node_id: IsogenyGraphNodeId(0)
        }
    );
}

#[test]
fn conservative_refinement_matches_node_local_evidence_without_incident_edges() {
    let conservative = refinement_at(0, CandidateRefinementStrategy::Conservative);
    let node_local = refinement_at(0, CandidateRefinementStrategy::NodeLocalLevelsOnly);

    assert_eq!(conservative, node_local);
}

#[test]
fn conservative_refinement_uses_only_unambiguous_edge_constraints() {
    let refinement = refinement_at(1, CandidateRefinementStrategy::Conservative);

    assert_only_node_or_unambiguous_edge_constraints(&refinement);
}

#[test]
fn conservative_refinement_runs_end_to_end_from_curve_graph() {
    let refinement = refinement_at(1, CandidateRefinementStrategy::Conservative);

    assert_refinement_candidates_stay_inside_initial_set(&refinement);
    assert!(refinement.constraints().iter().any(|constraint| matches!(
        constraint,
        LocalEndomorphismConstraint::NodeLevel {
            provenance: ConstraintSource::NodeReport {
                node_id: IsogenyGraphNodeId(0),
            },
            ..
        }
    )));
    assert_only_node_or_unambiguous_edge_constraints(&refinement);
}

#[test]
fn aggregate_refinement_report_refines_every_node_independently() {
    let (graph, report) = graph_report(1);

    let refinement_report = report
        .refine_candidates(CandidateRefinementStrategy::Conservative)
        .expect("aggregate refinement should build");

    assert_eq!(refinement_report.prime(), &BigUint::from(2u8));
    assert_eq!(
        refinement_report.strategy(),
        CandidateRefinementStrategy::Conservative
    );
    assert_eq!(
        refinement_report.node_refinements().len(),
        graph.node_count()
    );
    assert!(
        refinement_report
            .node_refinements()
            .iter()
            .all(|refinement| !refinement.initial_candidates().is_empty())
    );

    for refinement in refinement_report.node_refinements() {
        assert_eq!(
            refinement_report.refinement_for_node(refinement.node_id()),
            Some(refinement)
        );
        assert_refinement_candidates_stay_inside_initial_set(refinement);
        assert_only_node_or_unambiguous_edge_constraints(refinement);
    }

    assert!(
        refinement_report
            .refinement_for_node(IsogenyGraphNodeId(99))
            .is_none()
    );
}

#[test]
fn fixed_point_refinement_runs_end_to_end_from_curve_graph() {
    let (graph, report) = graph_report(1);

    let refinement_report = report
        .refine_candidates_to_fixed_point(CandidateRefinementStrategy::Conservative)
        .expect("fixed-point refinement should build");

    assert_eq!(refinement_report.prime(), &BigUint::from(2u8));
    assert_eq!(
        refinement_report.strategy(),
        CandidateRefinementStrategy::Conservative
    );
    assert_eq!(
        refinement_report.node_refinements().len(),
        graph.node_count()
    );
    for refinement in refinement_report.node_refinements() {
        assert_refinement_candidates_stay_inside_initial_set(refinement);
        assert_only_node_or_unambiguous_edge_constraints(refinement);
    }
}

#[test]
fn graph_report_records_observed_volcano_evidence_e2e() {
    let (_, report) = f17_graph_report(1, 2, 1);

    assert!(
        report
            .nodes()
            .iter()
            .any(|node| node.observed_allowed_levels().is_some())
    );
}

#[test]
fn fixed_point_refinement_can_reduce_candidates_e2e_from_a_real_curve() {
    let (_, report) = f17_graph_report(1, 2, 1);
    let fixed = report
        .fixed_point_candidate_refinement(CandidateRefinementStrategy::Conservative)
        .expect("fixed-point refinement should build");

    assert_eq!(fixed.rounds_with_eliminations(), 1);
    assert_eq!(total_initial_candidates(fixed.report()), 4);
    assert_eq!(total_surviving_candidates(fixed.report()), 2);
    assert!(fixed.report().node_refinements().iter().any(|node| {
        !node.eliminated_candidates().is_empty() && !node.surviving_candidates().is_empty()
    }));
}

#[test]
fn fixed_point_refinement_can_need_two_elimination_rounds_e2e() {
    let (_, report) = f17_graph_report(1, 0, 1);
    let fixed = report
        .fixed_point_candidate_refinement(CandidateRefinementStrategy::Conservative)
        .expect("fixed-point refinement should build");

    assert_eq!(fixed.rounds_with_eliminations(), 2);
    assert_eq!(total_initial_candidates(fixed.report()), 6);
    assert_eq!(total_surviving_candidates(fixed.report()), 0);
    assert!(fixed.report().node_refinements().iter().all(|node| {
        node.surviving_candidates().is_empty() && !node.eliminated_candidates().is_empty()
    }));
}

#[test]
fn refinement_rejects_missing_nodes() {
    let (_, report) = graph_report(0);

    let error = report
        .refine_candidates_for_node(
            IsogenyGraphNodeId(99),
            CandidateRefinementStrategy::NodeLocalLevelsOnly,
        )
        .expect_err("missing node should be rejected");

    assert_eq!(
        error,
        CandidateRefinementError::NodeNotFound {
            node_id: IsogenyGraphNodeId(99)
        }
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn fixed_point_refinement_preserves_graph_shape_and_candidate_subsets(
        a in 0i64..17,
        b in 0i64..17,
        depth in 0usize..=2,
    ) {
        let Ok(curve) = Curve17::new(F17::from_i64(a), F17::from_i64(b)) else {
            return Ok(());
        };
        let Ok(graph) = IsogenyGraphBuilder::new(curve, 2).max_depth(depth).build() else {
            return Ok(());
        };
        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("endomorphism report should build for a graph that built");
        let independent = report
            .refine_candidates(CandidateRefinementStrategy::Conservative)
            .expect("independent refinement should build");
        let fixed = report
            .refine_candidates_to_fixed_point(CandidateRefinementStrategy::Conservative)
            .expect("fixed-point refinement should build");

        prop_assert_eq!(fixed.node_refinements().len(), graph.node_count());
        prop_assert!(total_surviving_candidates(&fixed) <= total_initial_candidates(&fixed));
        assert_fixed_point_is_subset_of(&fixed, &independent);
        for refinement in fixed.node_refinements() {
            assert_refinement_candidates_stay_inside_initial_set(refinement);
        }
    }

    #[test]
    fn node_local_fixed_point_matches_independent_node_local_refinement(
        a in 0i64..17,
        b in 0i64..17,
        depth in 0usize..=2,
    ) {
        let Ok(curve) = Curve17::new(F17::from_i64(a), F17::from_i64(b)) else {
            return Ok(());
        };
        let Ok(graph) = IsogenyGraphBuilder::new(curve, 2).max_depth(depth).build() else {
            return Ok(());
        };
        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("endomorphism report should build for a graph that built");

        let independent = report
            .refine_candidates(CandidateRefinementStrategy::NodeLocalLevelsOnly)
            .expect("independent node-local refinement should build");
        let fixed = report
            .fixed_point_candidate_refinement(CandidateRefinementStrategy::NodeLocalLevelsOnly)
            .expect("fixed-point node-local refinement should build");

        prop_assert_eq!(fixed.rounds_with_eliminations(), 0);
        prop_assert_eq!(fixed.report(), &independent);
    }

    #[test]
    fn fixed_point_refinement_is_deterministic(
        a in 0i64..17,
        b in 0i64..17,
        depth in 0usize..=2,
    ) {
        let Ok(curve) = Curve17::new(F17::from_i64(a), F17::from_i64(b)) else {
            return Ok(());
        };
        let Ok(graph) = IsogenyGraphBuilder::new(curve, 2).max_depth(depth).build() else {
            return Ok(());
        };
        let report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("endomorphism report should build for a graph that built");

        let first = report
            .refine_candidates_to_fixed_point(CandidateRefinementStrategy::Conservative)
            .expect("first fixed-point refinement should build");
        let second = report
            .refine_candidates_to_fixed_point(CandidateRefinementStrategy::Conservative)
            .expect("second fixed-point refinement should build");

        prop_assert_eq!(first, second);
    }
}
