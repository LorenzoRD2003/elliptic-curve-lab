use num_bigint::BigUint;
use std::collections::{HashMap, HashSet};

use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::volcano_structure::{
        CraterShape, HorizontalEdgeReport, HorizontalEdgeStatus, VolcanoStructureReport,
    },
};

/// Crater report for a stored graph viewed as an ordinary `ℓ`-volcano.
///
/// The crater nodes are the certified surface `V₀` from the underlying
/// [`VolcanoStructureReport`]. Horizontal edges are certified by equal altitude
/// on that surface when possible, and weaker surface-like graph evidence is
/// kept in a separate status.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CraterReport {
    prime: BigUint,
    structure: VolcanoStructureReport,
    crater_nodes: Vec<IsogenyGraphNodeId>,
    horizontal_edges: Vec<HorizontalEdgeReport>,
    shape: CraterShape,
}

impl CraterReport {
    pub(crate) fn new(
        prime: BigUint,
        structure: VolcanoStructureReport,
        crater_nodes: Vec<IsogenyGraphNodeId>,
        horizontal_edges: Vec<HorizontalEdgeReport>,
        shape: CraterShape,
    ) -> Self {
        Self {
            prime,
            structure,
            crater_nodes,
            horizontal_edges,
            shape,
        }
    }

    /// Returns the chosen local prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the structural volcano report that supplies certified levels.
    pub fn structure(&self) -> &VolcanoStructureReport {
        &self.structure
    }

    /// Returns the certified crater nodes.
    pub fn nodes(&self) -> &[IsogenyGraphNodeId] {
        &self.crater_nodes
    }

    /// Returns edge reports that are horizontal, or possibly horizontal, in the
    /// crater evidence currently available.
    pub fn horizontal_edges(&self) -> &[HorizontalEdgeReport] {
        &self.horizontal_edges
    }

    /// Returns certified horizontal edges whose endpoints both lie in the crater.
    ///
    /// These are the graph edges available for deterministic crater walks:
    /// their status is [`HorizontalEdgeStatus::CertifiedByAltitude`], their
    /// source is a crater node, and their target is a crater node.
    pub fn certified_internal_horizontal_edges(&self) -> Vec<HorizontalEdgeReport> {
        let crater_nodes = self.crater_nodes.iter().copied().collect::<HashSet<_>>();

        self.horizontal_edges
            .iter()
            .filter(|edge| {
                edge.status() == HorizontalEdgeStatus::CertifiedByAltitude
                    && crater_nodes.contains(&edge.source())
                    && crater_nodes.contains(&edge.target())
            })
            .cloned()
            .collect()
    }

    /// Builds an outgoing-edge map keyed by the certified crater nodes.
    pub(crate) fn outgoing_edge_map(
        &self,
    ) -> HashMap<IsogenyGraphNodeId, Vec<HorizontalEdgeReport>> {
        self.crater_nodes
            .iter()
            .copied()
            .map(|node| (node, Vec::new()))
            .collect()
    }

    /// Returns the certified crater shape.
    pub fn shape(&self) -> CraterShape {
        self.shape
    }

    /// Returns the crater length when the certified shape determines it.
    pub fn crater_length(&self) -> Option<usize> {
        self.shape.crater_length()
    }

    /// Returns how many horizontal crater cycles are certified by the shape.
    pub fn horizontal_cycle_count(&self) -> usize {
        self.shape.horizontal_cycle_count()
    }

    /// Counts stored horizontal edges with the requested evidence status.
    pub fn horizontal_edge_count_by_status(&self, status: HorizontalEdgeStatus) -> usize {
        self.horizontal_edges
            .iter()
            .filter(|edge| edge.status() == status)
            .count()
    }
}
