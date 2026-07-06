use num_bigint::BigUint;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::endomorphisms::quadratic_orders::QuadraticDiscriminant;
use crate::isogenies::graphs::{
    EndomorphismVolcanoReport, IsogenyGraphBuilder, IsogenyGraphNodeId, VolcanoHeuristicComparison,
    VolcanoLikeLayering,
};

type F41 = crate::fields::Fp41;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

#[test]
fn report_can_compare_candidate_levels_with_a_precomputed_heuristic() {
    let candidate_set = QuadraticDiscriminant::new(-16)
        .factorization()
        .expect("-16 should factor canonically")
        .endomorphism_ring_candidates()
        .expect("candidate orders should construct");
    let heuristic = VolcanoLikeLayering::new(
        vec![vec![IsogenyGraphNodeId(0)], vec![IsogenyGraphNodeId(1)]],
        Vec::new(),
    );

    let report = EndomorphismVolcanoReport::from_graph_heuristic(
        &candidate_set,
        &BigUint::from(2u8),
        heuristic,
    )
    .expect("report should build");

    assert_eq!(report.possible_levels(), &[0, 1]);
    assert_eq!(report.heuristic_level_count(), 2);
    assert_eq!(report.local_order_candidate_count(), 2);
    assert_eq!(
        report.comparison_with_graph_heuristic(),
        &VolcanoHeuristicComparison::CompatibleLevelCount
    );
}

#[test]
fn report_marks_empty_heuristics_as_unavailable() {
    let candidate_set = QuadraticDiscriminant::new(-16)
        .factorization()
        .expect("-16 should factor canonically")
        .endomorphism_ring_candidates()
        .expect("candidate orders should construct");

    let report = EndomorphismVolcanoReport::from_graph_heuristic(
        &candidate_set,
        &BigUint::from(2u8),
        VolcanoLikeLayering::new(Vec::new(), Vec::new()),
    )
    .expect("report should build");

    assert_eq!(
        report.comparison_with_graph_heuristic(),
        &VolcanoHeuristicComparison::HeuristicUnavailable
    );
}

#[test]
fn report_can_run_the_existing_graph_heuristic_itself() {
    let candidate_set = QuadraticDiscriminant::new(-16)
        .factorization()
        .expect("-16 should factor canonically")
        .endomorphism_ring_candidates()
        .expect("candidate orders should construct");
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("depth-one graph should build");

    let report = candidate_set
        .volcano_report_from_graph_and_root(&BigUint::from(2u8), &graph, IsogenyGraphNodeId(0))
        .expect("report should build");

    assert_eq!(report.graph_heuristic().level_count(), 2);
    assert_eq!(
        report.comparison_with_graph_heuristic(),
        &VolcanoHeuristicComparison::CompatibleLevelCount
    );
}
