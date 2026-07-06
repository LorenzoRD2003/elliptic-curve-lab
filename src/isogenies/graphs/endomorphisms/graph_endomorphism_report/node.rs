use std::collections::BTreeSet;

use crate::elliptic_curves::endomorphisms::candidate_sets::{
    EndomorphismRingCandidateSet, VolcanoEndomorphismLevelCandidate,
};
use crate::isogenies::graphs::IsogenyGraphNodeId;

/// Endomorphism-side report for one stored graph node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IsogenyGraphEndomorphismNodeReport {
    node_id: IsogenyGraphNodeId,
    candidate_set: EndomorphismRingCandidateSet,
    local_levels: Vec<VolcanoEndomorphismLevelCandidate>,
    observed_allowed_levels: Option<BTreeSet<u32>>,
}

impl IsogenyGraphEndomorphismNodeReport {
    pub(crate) fn new(
        node_id: IsogenyGraphNodeId,
        candidate_set: EndomorphismRingCandidateSet,
        local_levels: Vec<VolcanoEndomorphismLevelCandidate>,
        observed_allowed_levels: Option<BTreeSet<u32>>,
    ) -> Self {
        Self {
            node_id,
            candidate_set,
            local_levels,
            observed_allowed_levels,
        }
    }

    /// Returns the node identifier.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns the Frobenius-compatible candidate orders for this node.
    pub fn candidate_set(&self) -> &EndomorphismRingCandidateSet {
        &self.candidate_set
    }

    /// Returns how many arithmetic `ℓ`-local level candidates were recorded
    /// for this node.
    pub fn local_level_candidate_count(&self) -> usize {
        self.local_levels.len()
    }

    /// Returns the distinct possible local levels for this node.
    pub fn possible_levels(&self) -> Vec<u32> {
        VolcanoEndomorphismLevelCandidate::distinct_levels_from(&self.local_levels)
    }

    /// Returns graph-observed local levels when the report has a conservative
    /// endpoint-role constraint for this node.
    ///
    /// These are still heuristic `ℓ`-volcano observations, not a certificate of
    /// the exact endomorphism ring. Absence means the node keeps the full
    /// arithmetic level set from `C₀`.
    pub fn observed_allowed_levels(&self) -> Option<&BTreeSet<u32>> {
        self.observed_allowed_levels.as_ref()
    }

    pub(crate) fn refinement_allowed_levels(&self) -> BTreeSet<u32> {
        self.observed_allowed_levels
            .clone()
            .unwrap_or_else(|| self.possible_levels().into_iter().collect())
    }
}
