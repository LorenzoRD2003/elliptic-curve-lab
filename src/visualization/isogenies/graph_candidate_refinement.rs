use crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder;
use crate::isogenies::graphs::endomorphisms::refinement::{
    CandidateEliminationReason, CandidateRefinementEdgeDirection, EndomorphismCandidateRefinement,
    IsogenyGraphCandidateRefinementReport, RefinementConfidence,
};
use crate::visualization::{
    Visualizable,
    shared::{comma_list, format_imaginary_quadratic_order_label as format_order_label},
};

/// Explains a graph-level candidate-refinement report in plain text.
///
/// This summarizes which orders `O_f` survived the observed local evidence at
/// the chosen prime `ℓ`. A unique survivor means “uniquely compatible with the
/// recorded evidence”, not a certificate that the order equals `End(E)`.
fn explain_graph_candidate_refinement_report(
    report: &IsogenyGraphCandidateRefinementReport,
) -> String {
    let total_initial = report
        .node_refinements()
        .iter()
        .map(|refinement| refinement.initial_candidates().len())
        .sum::<usize>();
    let total_surviving = report
        .node_refinements()
        .iter()
        .map(|refinement| refinement.surviving_candidates().len())
        .sum::<usize>();
    let total_eliminated = report
        .node_refinements()
        .iter()
        .map(|refinement| refinement.eliminated_candidates().len())
        .sum::<usize>();

    let mut lines = vec![
        "Endomorphism candidate refinement from graph evidence".to_string(),
        "-----------------------------------------------------".to_string(),
        format!("prime ℓ: {}", report.prime()),
        format!("strategy: {:?}", report.strategy()),
        format_refinement_confidence_line(report),
        format!("node refinements: {}", report.node_refinements().len()),
        format!("initial candidates total: {total_initial}"),
        format!("surviving candidates total: {total_surviving}"),
        format!("eliminated candidates total: {total_eliminated}"),
        "This report records compatibility with observed evidence; it does not certify exact End(E).".to_string(),
        String::new(),
        "Nodes:".to_string(),
    ];

    lines.extend(
        report
            .node_refinements()
            .iter()
            .flat_map(format_refinement_node_lines),
    );

    lines.join("\n")
}

impl Visualizable for IsogenyGraphCandidateRefinementReport {
    fn format_compact(&self) -> String {
        format!(
            "candidate refinement at ℓ = {} across {} nodes",
            self.prime(),
            self.node_refinements().len()
        )
    }

    fn describe(&self) -> String {
        explain_graph_candidate_refinement_report(self)
    }
}

fn format_refinement_node_lines(refinement: &EndomorphismCandidateRefinement) -> Vec<String> {
    let mut lines = vec![
        format!(
            "  v{}: initial {}, surviving {}, eliminated {}, constraints {}",
            refinement.node_id().0,
            refinement.initial_candidates().len(),
            refinement.surviving_candidates().len(),
            refinement.eliminated_candidates().len(),
            refinement.constraints().len(),
        ),
        format!(
            "    initial orders: {}",
            comma_list(
                refinement
                    .initial_candidates()
                    .candidate_orders()
                    .iter()
                    .map(format_order_label),
            )
        ),
        format!(
            "    surviving orders: {}",
            format_order_list(refinement.surviving_candidates())
        ),
    ];

    if let Some(unique_survivor) = refinement.unique_survivor() {
        lines.push(format!(
            "    unique survivor compatible with evidence: {}",
            format_order_label(unique_survivor)
        ));
    }

    if refinement.eliminated_candidates().is_empty() {
        lines.push("    eliminated orders: none".to_string());
    } else {
        lines.push("    eliminated orders:".to_string());
        lines.extend(
            refinement
                .eliminated_candidates()
                .iter()
                .map(|elimination| {
                    format!(
                        "      {}: {}",
                        format_order_label(elimination.candidate()),
                        format_elimination_reason(elimination.reason())
                    )
                }),
        );
    }

    lines
}

fn format_refinement_confidence_line(report: &IsogenyGraphCandidateRefinementReport) -> String {
    let confidence_values = distinct_refinement_confidences(report);
    match confidence_values.as_slice() {
        [] => "confidence: none".to_string(),
        [confidence] => format!("confidence: {confidence:?}"),
        many => format!(
            "confidence values: {}",
            comma_list(many.iter().map(|confidence| format!("{confidence:?}")))
        ),
    }
}

fn distinct_refinement_confidences(
    report: &IsogenyGraphCandidateRefinementReport,
) -> Vec<RefinementConfidence> {
    let mut confidence_values = Vec::new();
    for confidence in report
        .node_refinements()
        .iter()
        .map(|refinement| refinement.confidence())
    {
        if !confidence_values.contains(&confidence) {
            confidence_values.push(confidence);
        }
    }
    confidence_values
}

fn format_order_list(orders: &[ImaginaryQuadraticOrder]) -> String {
    if orders.is_empty() {
        "none".to_string()
    } else {
        comma_list(orders.iter().map(format_order_label))
    }
}

fn format_elimination_reason(reason: &CandidateEliminationReason) -> String {
    match reason {
        CandidateEliminationReason::IncompatibleLocalLevel {
            ell,
            candidate_level,
            allowed_levels,
        } => format!(
            "local level v_{}(f) = {} not in {:?}",
            ell, candidate_level, allowed_levels
        ),
        CandidateEliminationReason::IncompatibleIncidentEdgeRelation {
            ell,
            direction,
            adjacent_node,
            candidate_level,
            compatible_adjacent_levels,
            expected_relation,
        } => format!(
            "{} edge at ℓ = {} with v{}: local level {} incompatible with adjacent levels {:?} for {:?}",
            format_edge_direction(*direction),
            ell,
            adjacent_node.0,
            candidate_level,
            compatible_adjacent_levels,
            expected_relation
        ),
    }
}

fn format_edge_direction(direction: CandidateRefinementEdgeDirection) -> &'static str {
    match direction {
        CandidateRefinementEdgeDirection::Outgoing => "outgoing",
        CandidateRefinementEdgeDirection::Incoming => "incoming",
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::*;
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::isogenies::graphs::{
        IsogenyGraphBuilder, endomorphisms::refinement::CandidateRefinementStrategy,
    };
    use crate::visualization::traits::Visualizable;

    type F41 = crate::fields::Fp41;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn f41_curve() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn graph_candidate_refinement_explanation_mentions_evidence_and_non_certification() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");
        let endomorphism_report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("endomorphism report should build");
        let refinement_report = endomorphism_report
            .refine_candidates_to_fixed_point(CandidateRefinementStrategy::Conservative)
            .expect("candidate refinement should build");

        let explanation = explain_graph_candidate_refinement_report(&refinement_report);

        assert!(explanation.contains("Endomorphism candidate refinement from graph evidence"));
        assert!(explanation.contains("prime ℓ: 2"));
        assert!(explanation.contains("strategy: Conservative"));
        assert!(explanation.contains("confidence:"));
        assert!(explanation.contains("compatibility with observed evidence"));
        assert!(explanation.contains("does not certify exact End(E)"));
        assert!(explanation.contains("initial orders:"));
        assert!(explanation.contains("surviving orders:"));
    }

    #[test]
    fn graph_candidate_refinement_report_is_visualizable() {
        let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
            .max_depth(1)
            .build()
            .expect("depth-one graph should build");
        let endomorphism_report = graph
            .endomorphism_report_at(&BigUint::from(2u8))
            .expect("endomorphism report should build");
        let refinement_report = endomorphism_report
            .refine_candidates(CandidateRefinementStrategy::NodeLocalLevelsOnly)
            .expect("candidate refinement should build");

        assert_eq!(
            refinement_report.format_compact(),
            "candidate refinement at ℓ = 2 across 2 nodes"
        );
        assert!(
            refinement_report
                .describe()
                .contains("strategy: NodeLocalLevelsOnly")
        );
    }
}
