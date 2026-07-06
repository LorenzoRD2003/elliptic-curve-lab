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
