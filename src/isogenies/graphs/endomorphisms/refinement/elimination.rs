use std::collections::BTreeSet;

use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::quadratic_orders::ImaginaryQuadraticOrder;
use crate::isogenies::graphs::{
    IsogenyGraphNodeId, endomorphisms::IsogenyEdgeEndomorphismTentativeRelation,
};

/// Direction of an incident edge relative to the node being refined.
///
/// This keeps edge-based elimination reasons compact while preserving whether
/// the refined candidate was evaluated on the source side or target side of
/// the observed edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CandidateRefinementEdgeDirection {
    /// The refined node is the source of the edge.
    Outgoing,
    /// The refined node is the target of the edge.
    Incoming,
}

/// Structured reason why one candidate order did not survive a refinement run.
///
/// Reasons are deliberately data-shaped. This prevents callers from treating
/// vague prose as mathematical evidence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CandidateEliminationReason {
    /// The candidate's local level `v_ℓ(f)` is not one of the allowed node
    /// levels.
    IncompatibleLocalLevel {
        /// The chosen local prime `ℓ`.
        ell: BigUint,
        /// The local level of the eliminated candidate.
        candidate_level: u32,
        /// The set of levels compatible with the observed node evidence.
        allowed_levels: BTreeSet<u32>,
    },
    /// The candidate failed one incident edge relation.
    IncompatibleIncidentEdgeRelation {
        /// The chosen local prime `ℓ`.
        ell: BigUint,
        /// Whether the refined node is the source or target of the edge.
        direction: CandidateRefinementEdgeDirection,
        /// The other endpoint of the incident edge.
        adjacent_node: IsogenyGraphNodeId,
        /// The local level of the eliminated candidate on the refined node.
        candidate_level: u32,
        /// Local levels still compatible with the adjacent node.
        compatible_adjacent_levels: BTreeSet<u32>,
        /// The tentative relation expected along the edge.
        expected_relation: IsogenyEdgeEndomorphismTentativeRelation,
    },
}

/// One eliminated candidate order together with the first recorded reason.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CandidateElimination {
    candidate: ImaginaryQuadraticOrder,
    reason: CandidateEliminationReason,
}

impl CandidateElimination {
    pub(crate) fn new(
        candidate: ImaginaryQuadraticOrder,
        reason: CandidateEliminationReason,
    ) -> Self {
        Self { candidate, reason }
    }

    /// Returns the eliminated candidate order `O_f`.
    pub fn candidate(&self) -> &ImaginaryQuadraticOrder {
        &self.candidate
    }

    /// Returns the structured reason recorded for the elimination.
    pub fn reason(&self) -> &CandidateEliminationReason {
        &self.reason
    }
}
