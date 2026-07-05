use std::collections::BTreeSet;

use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::candidate_sets::EndomorphismRingCandidateSet;
use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{
        IsogenyGraphEndomorphismNodeReport,
        refinement::{
            CandidateElimination, CandidateEliminationReason, CandidateRefinementError,
            ConstraintSource, EndomorphismCandidateRefinement, IncidentEdgeRefinementConstraint,
            LocalEndomorphismConstraint, RefinementConfidence,
        },
    },
};

impl EndomorphismCandidateRefinement {
    pub(crate) fn from_node_local_levels(
        node_report: &IsogenyGraphEndomorphismNodeReport,
        ell: &BigUint,
    ) -> Result<Self, CandidateRefinementError> {
        let node_id = node_report.node_id();
        Self::from_node_level_constraint(
            node_id,
            node_report.candidate_set().clone(),
            ell.clone(),
            node_report.refinement_allowed_levels(),
            ConstraintSource::NodeReport { node_id },
        )
    }

    pub(crate) fn from_incident_unambiguous_edges(
        node_report: &IsogenyGraphEndomorphismNodeReport,
        ell: &BigUint,
        incident_constraints: Vec<IncidentEdgeRefinementConstraint>,
    ) -> Result<Self, CandidateRefinementError> {
        let node_id = node_report.node_id();
        let allowed_levels = node_report.refinement_allowed_levels();
        Self::from_constraints(
            node_id,
            node_report.candidate_set().clone(),
            ell.clone(),
            allowed_levels,
            ConstraintSource::NodeReport { node_id },
            incident_constraints,
        )
    }

    pub(crate) fn from_node_level_constraint(
        node_id: IsogenyGraphNodeId,
        initial_candidates: EndomorphismRingCandidateSet,
        ell: BigUint,
        allowed_levels: BTreeSet<u32>,
        provenance: ConstraintSource,
    ) -> Result<Self, CandidateRefinementError> {
        Self::from_constraints(
            node_id,
            initial_candidates,
            ell,
            allowed_levels,
            provenance,
            Vec::new(),
        )
    }

    pub(crate) fn from_constraints(
        node_id: IsogenyGraphNodeId,
        initial_candidates: EndomorphismRingCandidateSet,
        ell: BigUint,
        allowed_levels: BTreeSet<u32>,
        node_provenance: ConstraintSource,
        incident_constraints: Vec<IncidentEdgeRefinementConstraint>,
    ) -> Result<Self, CandidateRefinementError> {
        let mut surviving_candidates = Vec::new();
        let mut eliminated_candidates = Vec::new();

        for candidate in initial_candidates.candidate_orders() {
            let candidate_level = candidate
                .volcanic_level_at(&ell)
                .map_err(|_| CandidateRefinementError::InvalidLocalPrime)?
                .level();

            if allowed_levels.contains(&candidate_level) {
                if let Some(failing_constraint) = incident_constraints
                    .iter()
                    .find(|constraint| !constraint.allows_candidate_level(candidate_level))
                {
                    eliminated_candidates.push(CandidateElimination::new(
                        candidate.clone(),
                        failing_constraint.elimination_reason(&ell, candidate_level),
                    ));
                } else {
                    surviving_candidates.push(candidate.clone());
                }
            } else {
                eliminated_candidates.push(CandidateElimination::new(
                    candidate.clone(),
                    CandidateEliminationReason::IncompatibleLocalLevel {
                        ell: ell.clone(),
                        candidate_level,
                        allowed_levels: allowed_levels.clone(),
                    },
                ));
            }
        }

        let mut constraints = vec![LocalEndomorphismConstraint::NodeLevel {
            ell: ell.clone(),
            allowed_levels,
            provenance: node_provenance,
        }];
        constraints.extend(
            incident_constraints
                .iter()
                .map(|constraint| constraint.local_constraint(&ell)),
        );

        let confidence = if incident_constraints.is_empty() {
            RefinementConfidence::ConservativeLocalEvidence
        } else {
            RefinementConfidence::TentativeGraphEvidence
        };

        Ok(Self::new(
            node_id,
            initial_candidates,
            surviving_candidates,
            eliminated_candidates,
            constraints,
            confidence,
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use num_bigint::BigUint;

    use crate::elliptic_curves::endomorphisms::{
        candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
    };
    use crate::isogenies::graphs::{
        IsogenyGraphEdgeId, IsogenyGraphNodeId,
        endomorphisms::{
            IsogenyEdgeEndomorphismTentativeRelation,
            refinement::{
                CandidateEliminationReason, CandidateRefinementEdgeDirection, ConstraintSource,
                EndomorphismCandidateRefinement, IncidentEdgeRefinementConstraint,
                LocalEndomorphismConstraint, RefinementConfidence,
            },
        },
    };

    fn candidate_set(discriminant: i64) -> EndomorphismRingCandidateSet {
        QuadraticDiscriminant::new(discriminant)
            .factorization()
            .expect("test discriminant should factor canonically")
            .endomorphism_ring_candidates()
            .expect("candidate orders should construct")
    }

    #[test]
    fn node_level_constraint_records_incompatible_local_levels() {
        let candidates = candidate_set(-16);
        let allowed_levels = BTreeSet::from([0]);

        let refinement = EndomorphismCandidateRefinement::from_node_level_constraint(
            IsogenyGraphNodeId(7),
            candidates,
            BigUint::from(2u8),
            allowed_levels.clone(),
            ConstraintSource::NodeReport {
                node_id: IsogenyGraphNodeId(7),
            },
        )
        .expect("valid prime should refine candidates");

        assert_eq!(refinement.node_id(), IsogenyGraphNodeId(7));
        assert_eq!(refinement.surviving_candidates().len(), 1);
        assert_eq!(refinement.eliminated_candidates().len(), 1);
        assert_eq!(
            refinement.confidence(),
            RefinementConfidence::ConservativeLocalEvidence
        );

        assert_eq!(
            refinement.eliminated_candidates()[0].reason(),
            &CandidateEliminationReason::IncompatibleLocalLevel {
                ell: BigUint::from(2u8),
                candidate_level: 1,
                allowed_levels: allowed_levels.clone(),
            }
        );
        assert_eq!(refinement.constraints().len(), 1);
        assert_eq!(
            refinement.constraints()[0],
            LocalEndomorphismConstraint::NodeLevel {
                ell: BigUint::from(2u8),
                allowed_levels,
                provenance: ConstraintSource::NodeReport {
                    node_id: IsogenyGraphNodeId(7),
                },
            }
        );
    }

    #[test]
    fn incident_edge_constraint_records_incompatible_local_levels() {
        let candidates = candidate_set(-16);
        let adjacent_levels = BTreeSet::from([0]);

        let refinement = EndomorphismCandidateRefinement::from_constraints(
            IsogenyGraphNodeId(7),
            candidates,
            BigUint::from(2u8),
            BTreeSet::from([0, 1]),
            ConstraintSource::NodeReport {
                node_id: IsogenyGraphNodeId(7),
            },
            vec![IncidentEdgeRefinementConstraint::new(
                IsogenyGraphEdgeId(3),
                CandidateRefinementEdgeDirection::Outgoing,
                IsogenyGraphNodeId(7),
                IsogenyGraphNodeId(8),
                IsogenyGraphNodeId(8),
                adjacent_levels.clone(),
                IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal,
            )],
        )
        .expect("valid prime should refine candidates");

        assert_eq!(refinement.surviving_candidates().len(), 1);
        assert_eq!(refinement.eliminated_candidates().len(), 1);
        assert_eq!(
            refinement.confidence(),
            RefinementConfidence::TentativeGraphEvidence
        );
        assert_eq!(
            refinement.eliminated_candidates()[0].reason(),
            &CandidateEliminationReason::IncompatibleIncidentEdgeRelation {
                ell: BigUint::from(2u8),
                direction: CandidateRefinementEdgeDirection::Outgoing,
                adjacent_node: IsogenyGraphNodeId(8),
                candidate_level: 1,
                compatible_adjacent_levels: adjacent_levels,
                expected_relation: IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal,
            }
        );
        assert!(refinement.constraints().iter().any(|constraint| {
            matches!(
                constraint,
                LocalEndomorphismConstraint::EdgeRelation {
                    ell,
                    source_node: IsogenyGraphNodeId(7),
                    target_node: IsogenyGraphNodeId(8),
                    relation: IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal,
                    provenance: ConstraintSource::EdgeReport {
                        edge_id: IsogenyGraphEdgeId(3),
                    },
                } if ell == &BigUint::from(2u8)
            )
        }));
    }
}
