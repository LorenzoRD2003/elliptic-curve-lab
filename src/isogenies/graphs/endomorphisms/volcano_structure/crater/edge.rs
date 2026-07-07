use std::collections::HashSet;

use crate::isogenies::graphs::{IsogenyGraphEdgeId, IsogenyGraphNodeId};

/// Evidence status for one edge considered by the crater report.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HorizontalEdgeStatus {
    /// Both endpoints are certified crater nodes by equal altimeter level.
    CertifiedByAltitude,
    /// Weak graph evidence places both endpoints on a surface-like shell.
    SuspectedByWeakSurfaceEvidence,
    /// At least one endpoint lacks certified altitude because the graph is partial.
    NotCertifiableBecausePartialGraph,
}

impl HorizontalEdgeStatus {
    pub(crate) fn from_evidence(
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
        crater_nodes: &HashSet<IsogenyGraphNodeId>,
        weak_surface_nodes: &HashSet<IsogenyGraphNodeId>,
        partial_nodes: &HashSet<IsogenyGraphNodeId>,
    ) -> Option<Self> {
        if crater_nodes.contains(&source) && crater_nodes.contains(&target) {
            return Some(Self::CertifiedByAltitude);
        }

        let touches_partial = partial_nodes.contains(&source) || partial_nodes.contains(&target);
        if touches_partial
            && (crater_nodes.contains(&source)
                || crater_nodes.contains(&target)
                || weak_surface_nodes.contains(&source)
                || weak_surface_nodes.contains(&target))
        {
            return Some(Self::NotCertifiableBecausePartialGraph);
        }

        if weak_surface_nodes.contains(&source) && weak_surface_nodes.contains(&target) {
            return Some(Self::SuspectedByWeakSurfaceEvidence);
        }

        None
    }
}

/// Report for one stored edge that is horizontal, or possibly horizontal, in
/// the observed crater structure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HorizontalEdgeReport {
    edge_id: IsogenyGraphEdgeId,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    status: HorizontalEdgeStatus,
}

impl HorizontalEdgeReport {
    pub(crate) fn new(
        edge_id: IsogenyGraphEdgeId,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
        status: HorizontalEdgeStatus,
    ) -> Self {
        Self {
            edge_id,
            source,
            target,
            status,
        }
    }

    /// Returns the stored edge id.
    pub fn edge_id(&self) -> IsogenyGraphEdgeId {
        self.edge_id
    }

    /// Returns the source node of the stored directed edge.
    pub fn source(&self) -> IsogenyGraphNodeId {
        self.source
    }

    /// Returns the target node of the stored directed edge.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the evidence status for this crater-edge report.
    pub fn status(&self) -> HorizontalEdgeStatus {
        self.status
    }
}
