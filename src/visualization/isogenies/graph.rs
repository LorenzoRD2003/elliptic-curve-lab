use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphEdgeId, IsogenyGraphNodeId,
    IsogenyGraphVerificationReport, ReverseEdgeStatus, VolcanoLikeLayering, VolcanoRole,
    endomorphisms::IsogenyGraphEndomorphismReport,
};
use crate::visualization::shared::{comma_list, yes_no};
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

impl<C: GraphCurveModel> IsogenyGraph<C>
where
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
fn explain_isogeny_graph<C: GraphCurveModel + Visualizable>(graph: &IsogenyGraph<C>) -> String
where
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
        format!("has directed cycle: {}", yes_no(summary.has_directed_cycle)),
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
                        comma_list(level.iter().map(|node| format!("v{}", node.0)))
                    )
                }),
        );
    }

    lines.join("\n")
}

impl<C: GraphCurveModel + Visualizable> Visualizable for IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
    C::Elem: VisualizableField + Clone + Eq + Hash,
{
    fn format_compact(&self) -> String {
        let degree = self.edges().first().map(|edge| edge.degree()).unwrap_or(0);
        format!(
            "ℓ-isogeny graph: degree {}, {} nodes, {} edges",
            degree,
            self.node_count(),
            self.edge_count()
        )
    }

    fn describe(&self) -> String {
        explain_isogeny_graph(self)
    }
}

/// Explains one previously inferred volcano-like weak-BFS layering.
///
/// This helper reports the levels and node roles already present in `layers`;
/// it does not recompute or certify any arithmetic volcano structure.
fn explain_volcano_like_layers<C>(graph: &IsogenyGraph<C>, layers: &VolcanoLikeLayering) -> String
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
            comma_list(level.iter().map(|node| format!("v{}", node.0)))
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

impl Visualizable for VolcanoLikeLayering {
    fn format_compact(&self) -> String {
        format!(
            "volcano-like layering: {} levels, surface {}, middle {}, floor {}, isolated {}, unknown {}",
            self.level_count(),
            self.count_role(VolcanoRole::Surface),
            self.count_role(VolcanoRole::Middle),
            self.count_role(VolcanoRole::Floor),
            self.count_role(VolcanoRole::Isolated),
            self.count_role(VolcanoRole::Unknown),
        )
    }

    fn describe(&self) -> String {
        explain_stored_volcano_like_layering(self)
    }
}

/// Explains the exhaustive local verification report for a small isogeny graph.
///
/// The report is intentionally summarized instead of printed with `Debug`:
/// examples should show whether each verification family succeeded, and only
/// then provide compact reverse-edge details when they are informative.
fn explain_graph_verification_report(report: &IsogenyGraphVerificationReport) -> String {
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

impl Visualizable for IsogenyGraphVerificationReport {
    fn format_compact(&self) -> String {
        format!(
            "local graph verification: {}/{} maps land correctly, {}/{} reverse edges verified",
            self.edges_mapping_domain_to_codomain(),
            self.checked_edges(),
            self.reverse_edges_verified_as_dual(),
            self.checked_edges()
        )
    }

    fn describe(&self) -> String {
        explain_graph_verification_report(self)
    }
}

/// Explains the tentative endomorphism-side annotations attached to one graph.
///
/// This is a visualization of Frobenius-compatible candidate data only. It does
/// not certify the exact endomorphism ring of any node curve, and it does not
/// prove definitive horizontal/ascending/descending edge types.
fn explain_graph_endomorphism_report(report: &IsogenyGraphEndomorphismReport) -> String {
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

impl Visualizable for IsogenyGraphEndomorphismReport {
    fn format_compact(&self) -> String {
        format!(
            "endomorphism-side graph report at ℓ = {}: {} node reports, {} edge reports",
            self.prime(),
            self.nodes().len(),
            self.edges().len()
        )
    }

    fn describe(&self) -> String {
        explain_graph_endomorphism_report(self)
    }
}

/// Formats the graph as a compact adjacency list in dense node-id order.
fn format_adjacency_list<C: GraphCurveModel>(graph: &IsogenyGraph<C>) -> String
where
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
                format!("v{} -> {}", node.id().0, comma_list(targets))
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn summarize_volcano_heuristic<C: GraphCurveModel>(
    graph: &IsogenyGraph<C>,
) -> VolcanoHeuristicSummary
where
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

fn choose_volcano_root<C: GraphCurveModel>(graph: &IsogenyGraph<C>) -> Option<IsogenyGraphNodeId>
where
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

fn explain_stored_volcano_like_layering(layers: &VolcanoLikeLayering) -> String {
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
        "This stored layering records weak-BFS levels and node roles.".to_string(),
        String::new(),
        "Levels:".to_string(),
    ];

    if layers.levels().is_empty() {
        lines.push("  none".to_string());
    } else {
        lines.extend(layers.levels().iter().enumerate().map(|(index, level)| {
            format!(
                "  level {}: {}",
                index,
                comma_list(level.iter().map(|node| format!("v{}", node.0)))
            )
        }));
    }

    lines.push(String::new());
    lines.push("Node roles:".to_string());
    if layers.roles().is_empty() {
        lines.push("  none".to_string());
    } else {
        lines.extend(
            layers
                .roles()
                .iter()
                .map(|(node_id, role)| format!("  v{}: {:?}", node_id.0, role)),
        );
    }

    lines.join("\n")
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
mod tests;
