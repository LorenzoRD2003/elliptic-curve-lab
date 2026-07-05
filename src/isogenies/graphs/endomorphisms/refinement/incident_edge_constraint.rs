use std::collections::BTreeSet;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    IsogenyGraphEdgeId, IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyEdgeEndomorphismTentativeRelation, IsogenyGraphEndomorphismEdgeReport,
        refinement::{
            CandidateEliminationReason, CandidateRefinementEdgeDirection, ConstraintSource,
            LocalEndomorphismConstraint,
        },
    },
};

/// One unambiguous incident edge constraint used by a single-node refinement pass.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IncidentEdgeRefinementConstraint {
    edge_id: IsogenyGraphEdgeId,
    direction: CandidateRefinementEdgeDirection,
    source_node: IsogenyGraphNodeId,
    target_node: IsogenyGraphNodeId,
    adjacent_node: IsogenyGraphNodeId,
    adjacent_levels: BTreeSet<u32>,
    relation: IsogenyEdgeEndomorphismTentativeRelation,
}

impl IncidentEdgeRefinementConstraint {
    pub(crate) fn new(
        edge_id: IsogenyGraphEdgeId,
        direction: CandidateRefinementEdgeDirection,
        source_node: IsogenyGraphNodeId,
        target_node: IsogenyGraphNodeId,
        adjacent_node: IsogenyGraphNodeId,
        adjacent_levels: BTreeSet<u32>,
        relation: IsogenyEdgeEndomorphismTentativeRelation,
    ) -> Self {
        Self {
            edge_id,
            direction,
            source_node,
            target_node,
            adjacent_node,
            adjacent_levels,
            relation,
        }
    }

    pub(crate) fn from_edge_report(
        refined_node: IsogenyGraphNodeId,
        edge_report: &IsogenyGraphEndomorphismEdgeReport,
        adjacent_levels: BTreeSet<u32>,
    ) -> Option<Self> {
        let (direction, adjacent_node) = if edge_report.source() == refined_node {
            (
                CandidateRefinementEdgeDirection::Outgoing,
                edge_report.target(),
            )
        } else if edge_report.target() == refined_node {
            (
                CandidateRefinementEdgeDirection::Incoming,
                edge_report.source(),
            )
        } else {
            return None;
        };

        Some(Self::new(
            edge_report.edge_id(),
            direction,
            edge_report.source(),
            edge_report.target(),
            adjacent_node,
            adjacent_levels,
            edge_report.refinement_relation()?,
        ))
    }

    pub(crate) fn allows_candidate_level(&self, candidate_level: u32) -> bool {
        self.adjacent_levels.iter().any(|&adjacent_level| {
            let (source_level, target_level) = match self.direction {
                CandidateRefinementEdgeDirection::Outgoing => (candidate_level, adjacent_level),
                CandidateRefinementEdgeDirection::Incoming => (adjacent_level, candidate_level),
            };
            self.relation.allows_levels(source_level, target_level)
        })
    }

    pub(crate) fn local_constraint(&self, ell: &BigUint) -> LocalEndomorphismConstraint {
        LocalEndomorphismConstraint::EdgeRelation {
            ell: ell.clone(),
            source_node: self.source_node,
            target_node: self.target_node,
            relation: self.relation.clone(),
            provenance: ConstraintSource::EdgeReport {
                edge_id: self.edge_id,
            },
        }
    }

    pub(crate) fn elimination_reason(
        &self,
        ell: &BigUint,
        candidate_level: u32,
    ) -> CandidateEliminationReason {
        CandidateEliminationReason::IncompatibleIncidentEdgeRelation {
            ell: ell.clone(),
            direction: self.direction,
            adjacent_node: self.adjacent_node,
            candidate_level,
            compatible_adjacent_levels: self.adjacent_levels.clone(),
            expected_relation: self.relation.clone(),
        }
    }
}
