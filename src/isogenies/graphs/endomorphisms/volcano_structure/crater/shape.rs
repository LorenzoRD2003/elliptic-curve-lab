use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::volcano_structure::{HorizontalEdgeReport, HorizontalEdgeStatus},
};

/// Certified crater shape extracted from horizontal surface edges.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CraterShape {
    /// No certified crater nodes were found.
    EmptyCertifiedCrater,
    /// The certified crater has one node. Self-loops are counted as stored
    /// directed horizontal edges on that node.
    Singleton { self_loop_count: usize },
    /// The certified crater has two nodes. The count records stored directed
    /// horizontal edges between them.
    TwoVertex { directed_edge_count: usize },
    /// The certified crater is a simple cycle on at least three nodes.
    Cycle { length: usize },
    /// The certified crater nodes or horizontal edges do not yet determine one
    /// of the standard crater shapes.
    PartialOrAmbiguous,
}

impl CraterShape {
    pub(crate) fn from_crater_evidence(
        crater_nodes: &[IsogenyGraphNodeId],
        horizontal_edges: &[HorizontalEdgeReport],
    ) -> Self {
        let certified_edges = horizontal_edges
            .iter()
            .filter(|edge| edge.status() == HorizontalEdgeStatus::CertifiedByAltitude)
            .collect::<Vec<_>>();

        match crater_nodes.len() {
            0 => Self::EmptyCertifiedCrater,
            1 => {
                let node = crater_nodes[0];
                let self_loop_count = certified_edges
                    .iter()
                    .filter(|edge| edge.source() == node && edge.target() == node)
                    .count();
                Self::Singleton { self_loop_count }
            }
            2 => {
                let nodes = crater_nodes.iter().copied().collect::<HashSet<_>>();
                let directed_edge_count = certified_edges
                    .iter()
                    .filter(|edge| {
                        edge.source() != edge.target()
                            && nodes.contains(&edge.source())
                            && nodes.contains(&edge.target())
                    })
                    .count();
                Self::TwoVertex {
                    directed_edge_count,
                }
            }
            length if certified_edges_form_simple_cycle(crater_nodes, &certified_edges) => {
                Self::Cycle { length }
            }
            _ => Self::PartialOrAmbiguous,
        }
    }

    /// Returns the crater length when the certified shape determines it.
    pub fn crater_length(&self) -> Option<usize> {
        match self {
            Self::EmptyCertifiedCrater | Self::PartialOrAmbiguous => None,
            Self::Singleton { .. } => Some(1),
            Self::TwoVertex { .. } => Some(2),
            Self::Cycle { length } => Some(*length),
        }
    }

    /// Returns how many horizontal crater cycles this shape certifies.
    pub fn horizontal_cycle_count(&self) -> usize {
        match self {
            Self::Singleton { .. } | Self::TwoVertex { .. } | Self::Cycle { .. } => 1,
            Self::EmptyCertifiedCrater | Self::PartialOrAmbiguous => 0,
        }
    }
}

fn certified_edges_form_simple_cycle(
    crater_nodes: &[IsogenyGraphNodeId],
    certified_edges: &[&HorizontalEdgeReport],
) -> bool {
    let crater_node_set = crater_nodes.iter().copied().collect::<HashSet<_>>();
    let mut neighbors = crater_nodes
        .iter()
        .copied()
        .map(|node| (node, BTreeSet::new()))
        .collect::<HashMap<_, _>>();

    for edge in certified_edges {
        if edge.source() == edge.target() {
            return false;
        }
        if !crater_node_set.contains(&edge.source()) || !crater_node_set.contains(&edge.target()) {
            return false;
        }
        neighbors
            .get_mut(&edge.source())
            .expect("source is a crater node")
            .insert(edge.target());
        neighbors
            .get_mut(&edge.target())
            .expect("target is a crater node")
            .insert(edge.source());
    }

    if neighbors.values().any(|adjacent| adjacent.len() != 2) {
        return false;
    }

    connected_crater_node_count(crater_nodes[0], &neighbors) == crater_nodes.len()
}

fn connected_crater_node_count(
    start: IsogenyGraphNodeId,
    neighbors: &HashMap<IsogenyGraphNodeId, BTreeSet<IsogenyGraphNodeId>>,
) -> usize {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::from([start]);

    while let Some(node) = queue.pop_front() {
        if !visited.insert(node) {
            continue;
        }
        if let Some(adjacent) = neighbors.get(&node) {
            queue.extend(adjacent.iter().copied());
        }
    }
    visited.len()
}
