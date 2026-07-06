use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphEdgeId, IsogenyGraphNodeId,
    IsogenyGraphVerificationReport, ReverseEdgeStatus, VolcanoLikeLayering, VolcanoRole,
    endomorphisms::IsogenyGraphEndomorphismReport,
};
use crate::visualization::{Visualizable, VisualizableField};

/// Root-dependent educational volcano heuristic attached to one graph summary.
///
/// The root is chosen deterministically as the smallest dense node id inside
/// the largest weakly connected component. The result is stable for one stored
/// graph, but remains only a graph-theoretic teaching aid.
#[derive(Clone, Debug)]
pub struct VolcanoHeuristicSummary {
    pub root: Option<IsogenyGraphNodeId>,
    pub levels: Vec<Vec<IsogenyGraphNodeId>>,
    pub surface_nodes: usize,
    pub middle_nodes: usize,
    pub floor_nodes: usize,
    pub isolated_nodes: usize,
    pub unknown_nodes: usize,
}

impl VolcanoHeuristicSummary {
    pub fn levels(&self) -> &[Vec<IsogenyGraphNodeId>] {
        &self.levels
    }

    pub fn level_count(&self) -> usize {
        self.levels.len()
    }

    pub fn is_empty(&self) -> bool {
        self.levels.is_empty()
    }
}

/// Compact structural summary of one educational `ℓ`-isogeny graph.
#[derive(Clone, Debug)]
pub struct IsogenyGraphSummary {
    pub node_count: usize,
    pub edge_count: usize,
    /// Current implementation simplification: this is read from the first stored
    /// edge and therefore assumes the whole graph was built for one fixed
    /// prime degree `ℓ`.
    ///
    /// Technical debt:
    /// if the graph container later allows mixed degrees, this field will need
    /// to become a richer summary surface such as a degree set, a histogram,
    /// or an optional “uniform degree” marker instead of one plain `usize`.
    pub degree: usize,
    pub connected_component_count: usize,
    pub has_directed_cycle: bool,
    pub self_loops: usize,
    pub repeated_j_invariants: usize,
    pub min_out_degree: usize,
    pub max_out_degree: usize,
    /// Educational, root-dependent weak-BFS layering inspired by volcano-like
    /// pictures of small `ℓ`-isogeny graphs.
    ///
    /// Technical debt:
    /// this is still only a graph-only heuristic and does not certify a true
    /// arithmetic isogeny volcano.
    pub volcano_like: VolcanoHeuristicSummary,
}

impl<C> IsogenyGraph<C>
where
    C: GraphCurveModel,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
    C::Elem: Clone + Eq + Hash,
{
    /// Returns a small structural summary of the stored graph.
    ///
    /// The component count is computed on the underlying undirected graph,
    /// since the educational question is usually “how many pieces does the
    /// explored graph break into?” rather than strong connectivity.
    ///
    /// The `degree` field assumes a same-`ℓ` graph. When no edges are present
    /// it is reported as `0`.
    ///
    /// Technical debt:
    /// this summary currently does not try to detect or explain mixed-degree
    /// edge sets inside one graph container.
    ///
    /// The `volcano_like` field is intentionally heuristic and root-dependent.
    /// It summarizes one weak-BFS layering chosen from the largest weakly
    /// connected component, not a proof of a true arithmetic volcano.
    pub fn summary(&self) -> IsogenyGraphSummary {
        let node_count = self.node_count();
        let edge_count = self.edge_count();
        let degree = self.edges().first().map(|edge| edge.degree()).unwrap_or(0);
        let self_loops = self
            .edges()
            .iter()
            .filter(|edge| edge.source() == edge.target())
            .count();
        let unique_j_invariants = self
            .nodes()
            .iter()
            .map(|node| node.j_invariant())
            .collect::<HashSet<_>>()
            .len();
        let repeated_j_invariants = node_count.saturating_sub(unique_j_invariants);

        let (min_out_degree, max_out_degree) = if self.nodes().is_empty() {
            (0, 0)
        } else {
            let mut out_degrees = self.nodes().iter().map(|node| self.out_degree(node.id()));
            let first = out_degrees
                .next()
                .expect("non-empty node set yields at least one out degree");
            out_degrees.fold((first, first), |(min_degree, max_degree), degree| {
                (min_degree.min(degree), max_degree.max(degree))
            })
        };
        let volcano_like = summarize_volcano_heuristic(self);

        IsogenyGraphSummary {
            node_count,
            edge_count,
            degree,
            connected_component_count: self.weakly_connected_components().len(),
            has_directed_cycle: self.has_directed_cycle(),
            self_loops,
            repeated_j_invariants,
            min_out_degree,
            max_out_degree,
            volcano_like,
        }
    }
}

/// Explains an educational `ℓ`-isogeny graph in plain text.
pub fn explain_isogeny_graph<C>(graph: &IsogenyGraph<C>) -> String
where
    C: GraphCurveModel + Visualizable,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
    C::Elem: VisualizableField + Clone + Eq + Hash,
{
    let summary = graph.summary();
    let mut lines = vec![
        "ℓ-isogeny graph summary".to_string(),
        "-----------------------".to_string(),
        format!("degree ℓ: {}", summary.degree),
        format!("nodes: {}", summary.node_count),
        format!("edges: {}", summary.edge_count),
        format!(
            "connected components: {}",
            summary.connected_component_count
        ),
        format!(
            "has directed cycle: {}",
            if summary.has_directed_cycle {
                "yes"
            } else {
                "no"
            }
        ),
        format!("self loops: {}", summary.self_loops),
        format!("repeated j-invariants: {}", summary.repeated_j_invariants),
        format!(
            "volcano-like root: {}",
            summary
                .volcano_like
                .root
                .map(|root| format!("v{}", root.0))
                .unwrap_or_else(|| "none".to_string())
        ),
        format!(
            "volcano-like levels: {}",
            summary.volcano_like.level_count()
        ),
        format!(
            "volcano-like roles: surface {}, middle {}, floor {}, isolated {}, unknown {}",
            summary.volcano_like.surface_nodes,
            summary.volcano_like.middle_nodes,
            summary.volcano_like.floor_nodes,
            summary.volcano_like.isolated_nodes,
            summary.volcano_like.unknown_nodes,
        ),
        String::new(),
        "Nodes:".to_string(),
    ];

    lines.extend(graph.nodes().iter().map(|node| {
        format!(
            "  v{}: j = {}, curve = {}",
            node.id().0,
            node.j_invariant().format_elem(),
            node.representative().format_compact()
        )
    }));

    lines.push(String::new());
    lines.push("Edges:".to_string());
    lines.extend(graph.edges().iter().map(|edge| {
        format!(
            "  e{}: v{} -> v{}, degree {}, kernel size {}",
            edge.id().0,
            edge.source().0,
            edge.target().0,
            edge.degree(),
            edge.kernel_order()
        )
    }));

    lines.push(String::new());
    lines.push("Adjacency list:".to_string());
    lines.push(format_adjacency_list(graph));

    if !summary.volcano_like.is_empty() {
        lines.push(String::new());
        lines.push("Volcano-like levels (heuristic):".to_string());
        lines.extend(
            summary
                .volcano_like
                .levels()
                .iter()
                .enumerate()
                .map(|(index, level)| {
                    format!(
                        "  level {}: {}",
                        index,
                        level
                            .iter()
                            .map(|node| format!("v{}", node.0))
                            .collect::<Vec<_>>()
                            .join(", ")
                    )
                }),
        );
    }

    lines.join("\n")
}

/// Explains one previously inferred volcano-like weak-BFS layering.
///
/// This helper reports the levels and node roles already present in `layers`;
/// it does not recompute or certify any arithmetic volcano structure.
pub fn explain_volcano_like_layers<C>(
    graph: &IsogenyGraph<C>,
    layers: &VolcanoLikeLayering,
) -> String
where
    C: GraphCurveModel + Visualizable,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
    C::Elem: VisualizableField + Clone + Eq + Hash,
{
    let mut lines = vec![
        "Volcano-like layering (heuristic)".to_string(),
        "--------------------------------".to_string(),
        format!("levels: {}", layers.level_count()),
        format!(
            "roles: surface {}, middle {}, floor {}, isolated {}, unknown {}",
            layers.count_role(VolcanoRole::Surface),
            layers.count_role(VolcanoRole::Middle),
            layers.count_role(VolcanoRole::Floor),
            layers.count_role(VolcanoRole::Isolated),
            layers.count_role(VolcanoRole::Unknown),
        ),
        "This layering is weak-BFS-based and educational only.".to_string(),
        String::new(),
        "Levels:".to_string(),
    ];

    lines.extend(layers.levels().iter().enumerate().map(|(index, level)| {
        format!(
            "  level {}: {}",
            index,
            level
                .iter()
                .map(|node| format!("v{}", node.0))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }));

    lines.push(String::new());
    lines.push("Node roles:".to_string());
    lines.extend(layers.roles().iter().map(|(node_id, role)| {
        let node = graph
            .node(*node_id)
            .expect("layering should only reference existing nodes");
        format!(
            "  v{}: {:?}, j = {}, curve = {}",
            node_id.0,
            role,
            node.j_invariant().format_elem(),
            node.representative().format_compact()
        )
    }));

    lines.join("\n")
}

/// Explains the exhaustive local verification report for a small isogeny graph.
///
/// The report is intentionally summarized instead of printed with `Debug`:
/// examples should show whether each verification family succeeded, and only
/// then provide compact reverse-edge details when they are informative.
pub fn explain_graph_verification_report(report: &IsogenyGraphVerificationReport) -> String {
    let reverse_status_counts = count_reverse_edge_statuses(report.reverse_edge_statuses());
    let mut lines = vec![
        "Local graph verification report".to_string(),
        "-------------------------------".to_string(),
        format!("checked edges: {}", report.checked_edges()),
        format!(
            "maps domain to codomain: {}/{}",
            report.edges_mapping_domain_to_codomain(),
            report.checked_edges()
        ),
        format!(
            "maps kernel to identity: {}/{}",
            report.edges_mapping_kernel_to_identity(),
            report.checked_edges()
        ),
        format!(
            "homomorphism law verified: {}/{}",
            report.edges_homomorphism_ok(),
            report.checked_edges()
        ),
        format!(
            "reverse edges verified as dual: {}/{}",
            report.reverse_edges_verified_as_dual(),
            report.checked_edges()
        ),
        format!(
            "reverse-edge statuses: verified {}, present-not-verified {}, missing {}",
            reverse_status_counts.verified_as_dual,
            reverse_status_counts.present_but_not_verified,
            reverse_status_counts.missing
        ),
    ];

    if reverse_edge_details_are_informative(report.reverse_edge_statuses()) {
        lines.push(String::new());
        lines.push("Reverse-edge details:".to_string());
        lines.extend(
            report
                .reverse_edge_statuses()
                .iter()
                .map(|(edge_id, status)| {
                    format!("  e{}: {}", edge_id.0, format_reverse_edge_status(*status))
                }),
        );
    }

    lines.join("\n")
}

/// Explains the tentative endomorphism-side annotations attached to one graph.
///
/// This is a visualization of Frobenius-compatible candidate data only. It does
/// not certify the exact endomorphism ring of any node curve, and it does not
/// prove definitive horizontal/ascending/descending edge types.
pub fn explain_graph_endomorphism_report(report: &IsogenyGraphEndomorphismReport) -> String {
    let mut lines = vec![
        "Tentative endomorphism-side report".to_string(),
        "----------------------------------".to_string(),
        format!("prime ℓ: {}", report.prime()),
        format!("node reports: {}", report.nodes().len()),
        format!("edge reports: {}", report.edges().len()),
        "This report is Frobenius-compatible only; it does not certify exact End(E).".to_string(),
        String::new(),
        "Nodes:".to_string(),
    ];

    lines.extend(report.nodes().iter().map(|node| {
        format!(
            "  v{}: candidates {}, possible levels {:?}",
            node.node_id().0,
            node.candidate_set().len(),
            node.possible_levels()
        )
    }));

    lines.push(String::new());
    lines.push("Edges:".to_string());
    lines.extend(report.edges().iter().map(|edge| {
        let relation = edge.relation();
        format!(
            "  e{}: v{} -> v{}, {:?}, source levels {:?}, target levels {:?}",
            edge.edge_id().0,
            edge.source().0,
            edge.target().0,
            relation.relation(),
            relation.source_possible_levels(),
            relation.target_possible_levels()
        )
    }));

    lines.join("\n")
}

/// Formats the graph as a compact adjacency list in dense node-id order.
pub fn format_adjacency_list<C>(graph: &IsogenyGraph<C>) -> String
where
    C: GraphCurveModel,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    graph
        .nodes()
        .iter()
        .map(|node| {
            let targets = graph
                .outgoing_edges(node.id())
                .map(|edge| format!("v{}", edge.target().0))
                .collect::<Vec<_>>();

            if targets.is_empty() {
                format!("v{} ->", node.id().0)
            } else {
                format!("v{} -> {}", node.id().0, targets.join(", "))
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn summarize_volcano_heuristic<C>(graph: &IsogenyGraph<C>) -> VolcanoHeuristicSummary
where
    C: GraphCurveModel,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    let Some(root) = choose_volcano_root(graph) else {
        return VolcanoHeuristicSummary {
            root: None,
            levels: Vec::new(),
            surface_nodes: 0,
            middle_nodes: 0,
            floor_nodes: 0,
            isolated_nodes: 0,
            unknown_nodes: 0,
        };
    };

    let layering = graph.infer_volcano_like_layers(root);
    let (surface_nodes, middle_nodes, floor_nodes, isolated_nodes, unknown_nodes) =
        count_volcano_roles(&layering);

    VolcanoHeuristicSummary {
        root: Some(root),
        levels: layering.levels().to_vec(),
        surface_nodes,
        middle_nodes,
        floor_nodes,
        isolated_nodes,
        unknown_nodes,
    }
}

fn choose_volcano_root<C>(graph: &IsogenyGraph<C>) -> Option<IsogenyGraphNodeId>
where
    C: GraphCurveModel,
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    graph
        .weakly_connected_components()
        .into_iter()
        .max_by(|left, right| left.len().cmp(&right.len()).then_with(|| right.cmp(left)))
        .and_then(|component| component.into_iter().min())
}

fn count_volcano_roles(layering: &VolcanoLikeLayering) -> (usize, usize, usize, usize, usize) {
    layering.roles().iter().fold(
        (0, 0, 0, 0, 0),
        |(surface, middle, floor, isolated, unknown), (_, role)| match role {
            VolcanoRole::Surface => (surface + 1, middle, floor, isolated, unknown),
            VolcanoRole::Middle => (surface, middle + 1, floor, isolated, unknown),
            VolcanoRole::Floor => (surface, middle, floor + 1, isolated, unknown),
            VolcanoRole::Isolated => (surface, middle, floor, isolated + 1, unknown),
            VolcanoRole::Unknown => (surface, middle, floor, isolated, unknown + 1),
        },
    )
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct ReverseEdgeStatusCounts {
    missing: usize,
    present_but_not_verified: usize,
    verified_as_dual: usize,
}

fn count_reverse_edge_statuses(
    statuses: &[(IsogenyGraphEdgeId, ReverseEdgeStatus)],
) -> ReverseEdgeStatusCounts {
    statuses.iter().fold(
        ReverseEdgeStatusCounts::default(),
        |mut counts, (_, status)| {
            match status {
                ReverseEdgeStatus::Missing => counts.missing += 1,
                ReverseEdgeStatus::PresentButNotVerifiedAsDual => {
                    counts.present_but_not_verified += 1;
                }
                ReverseEdgeStatus::VerifiedAsDual => counts.verified_as_dual += 1,
            }
            counts
        },
    )
}

fn reverse_edge_details_are_informative(
    statuses: &[(IsogenyGraphEdgeId, ReverseEdgeStatus)],
) -> bool {
    let Some((_, first_status)) = statuses.first() else {
        return false;
    };

    statuses.iter().any(|(_, status)| status != first_status)
}

fn format_reverse_edge_status(status: ReverseEdgeStatus) -> &'static str {
    match status {
        ReverseEdgeStatus::Missing => "missing",
        ReverseEdgeStatus::PresentButNotVerifiedAsDual => "present, not verified as dual",
        ReverseEdgeStatus::VerifiedAsDual => "verified as dual",
    }
}

#[cfg(test)]
mod tests {

    use std::collections::HashSet;

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::traits::Field;
    use crate::isogenies::graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId};
    use crate::visualization::isogenies::{
        IsogenyGraphSummary, VolcanoHeuristicSummary, explain_graph_endomorphism_report,
        explain_graph_verification_report, explain_isogeny_graph, explain_volcano_like_layers,
        format_adjacency_list,
    };
    use num_bigint::BigUint;

    type F5 = crate::fields::Fp5;
    type F41 = crate::fields::Fp41;
    type Curve41 = ShortWeierstrassCurve<F41>;
    type Curve5 = ShortWeierstrassCurve<F5>;

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
}
