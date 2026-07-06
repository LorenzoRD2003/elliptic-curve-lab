use crate::isogenies::graphs::{
    IsogenyGraphEdgeId, IsogenyGraphNodeId,
    endomorphisms::{IsogenyEdgeEndomorphismReport, IsogenyEdgeEndomorphismTentativeRelation},
};

/// Tentative endomorphism-side report for one stored graph edge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismEdgeReport {
    edge_id: IsogenyGraphEdgeId,
    source: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    relation: IsogenyEdgeEndomorphismReport,
    observed_relation: Option<IsogenyEdgeEndomorphismTentativeRelation>,
}

impl IsogenyGraphEndomorphismEdgeReport {
    pub(crate) fn new(
        edge_id: IsogenyGraphEdgeId,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
        relation: IsogenyEdgeEndomorphismReport,
        observed_relation: Option<IsogenyEdgeEndomorphismTentativeRelation>,
    ) -> Self {
        Self {
            edge_id,
            source,
            target,
            relation,
            observed_relation,
        }
    }

    /// Returns the edge identifier.
    pub fn edge_id(&self) -> IsogenyGraphEdgeId {
        self.edge_id
    }

    /// Returns the source node identifier.
    pub fn source(&self) -> IsogenyGraphNodeId {
        self.source
    }

    /// Returns the target node identifier.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the tentative endomorphism-side edge relation report.
    pub fn relation(&self) -> &IsogenyEdgeEndomorphismReport {
        &self.relation
    }

    /// Returns the graph-observed relation used by refinement when present.
    ///
    /// This relation comes from surface-anchored weak-BFS volcano evidence and
    /// is intentionally separate from the arithmetic candidate-set relation.
    pub fn observed_relation(&self) -> Option<&IsogenyEdgeEndomorphismTentativeRelation> {
        self.observed_relation.as_ref()
    }

    pub(crate) fn refinement_relation(&self) -> Option<IsogenyEdgeEndomorphismTentativeRelation> {
        self.observed_relation
            .clone()
            .or_else(|| self.relation.relation().as_unambiguous())
    }
}
