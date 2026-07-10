use std::collections::BTreeMap;

use super::*;
use crate::isogenies::{
    class_group_action::{
        CraterDirectionCertification, CraterOrientationWitness, CraterOrientationWitnessError,
    },
    graphs::{
        IsogenyGraphBuilder, IsogenyGraphEdgeId, IsogenyGraphNodeId,
        endomorphisms::{HorizontalEdgeReport, HorizontalEdgeStatus},
    },
};

fn edge(id: usize, source: usize, target: usize) -> HorizontalEdgeReport {
    HorizontalEdgeReport::new(
        IsogenyGraphEdgeId(id),
        IsogenyGraphNodeId(source),
        IsogenyGraphNodeId(target),
        HorizontalEdgeStatus::CertifiedByAltitude,
    )
}

fn two_node_crater() -> crate::isogenies::graphs::endomorphisms::CraterReport {
    crater_report_with_nodes(
        bu(3),
        vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)],
        vec![edge(0, 0, 1), edge(1, 1, 0)],
    )
}

fn orientation_01() -> BTreeMap<IsogenyGraphNodeId, IsogenyGraphNodeId> {
    BTreeMap::from([
        (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
        (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
    ])
}

#[test]
fn user_orientation_witness_accepts_a_certified_crater_cycle() {
    let crater = two_node_crater();

    let witness = CraterOrientationWitness::new(&crater, orientation_01())
        .expect("certified two-node crater cycle should orient");

    assert_eq!(
        witness.successor(IsogenyGraphNodeId(0)),
        Some(IsogenyGraphNodeId(1))
    );
    assert_eq!(
        witness.oriented_cycle_from(IsogenyGraphNodeId(0)),
        Some(vec![
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(0)
        ])
    );
}

#[test]
fn user_orientation_witness_rejects_missing_successor() {
    let crater = two_node_crater();
    let successors = BTreeMap::from([(IsogenyGraphNodeId(0), IsogenyGraphNodeId(1))]);

    let error = CraterOrientationWitness::new(&crater, successors)
        .expect_err("every crater node should have a declared successor");

    assert_eq!(
        error,
        CraterOrientationWitnessError::MissingSuccessor {
            source: IsogenyGraphNodeId(1)
        }
    );
}

#[test]
fn user_orientation_witness_rejects_source_outside_crater() {
    let crater = two_node_crater();
    let mut successors = orientation_01();
    successors.insert(IsogenyGraphNodeId(99), IsogenyGraphNodeId(0));

    let error = CraterOrientationWitness::new(&crater, successors)
        .expect_err("sources should be certified crater nodes");

    assert_eq!(
        error,
        CraterOrientationWitnessError::SourceOutsideCrater {
            source: IsogenyGraphNodeId(99)
        }
    );
}

#[test]
fn user_orientation_witness_rejects_target_outside_crater() {
    let crater = two_node_crater();
    let successors = BTreeMap::from([
        (IsogenyGraphNodeId(0), IsogenyGraphNodeId(99)),
        (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
    ]);

    let error = CraterOrientationWitness::new(&crater, successors)
        .expect_err("targets should be certified crater nodes");

    assert_eq!(
        error,
        CraterOrientationWitnessError::TargetOutsideCrater {
            source: IsogenyGraphNodeId(0),
            target: IsogenyGraphNodeId(99)
        }
    );
}

#[test]
fn user_orientation_witness_rejects_missing_certified_horizontal_edge() {
    let crater = crater_report_with_nodes(
        bu(3),
        vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)],
        Vec::new(),
    );

    let error = CraterOrientationWitness::new(&crater, orientation_01())
        .expect_err("successors should follow certified horizontal crater edges");

    assert_eq!(
        error,
        CraterOrientationWitnessError::MissingCertifiedHorizontalEdge {
            source: IsogenyGraphNodeId(0),
            target: IsogenyGraphNodeId(1)
        }
    );
}

#[test]
fn user_orientation_witness_inverse_reverses_the_successors() {
    let crater = two_node_crater();
    let witness = CraterOrientationWitness::new(&crater, orientation_01())
        .expect("certified two-node crater cycle should orient");

    let inverse = witness
        .inverse(&crater)
        .expect("inverse orientation should validate against the same crater");

    assert_eq!(
        inverse.oriented_cycle_from(IsogenyGraphNodeId(1)),
        Some(vec![
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1)
        ])
    );
}

#[test]
fn labeled_walk_can_be_wrapped_with_user_orientation() {
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");
    let ideal = split_three_ideal();
    let crater = graph
        .volcano_crater_report(ideal.norm())
        .expect("crater report should build for the ideal norm");
    let labeled = graph
        .labeled_crater_walk_report(&class_group_minus_23(), ideal, IsogenyGraphNodeId(0))
        .expect("labeled crater walk should build");
    let witness = CraterOrientationWitness::new(&crater, orientation_01())
        .expect("certified two-node crater cycle should orient");

    let oriented = labeled
        .with_user_orientation(witness)
        .expect("witness should attach to the labeled walk");

    assert_eq!(
        oriented.direction_certification(),
        CraterDirectionCertification::UserSuppliedArithmeticOrientation
    );
    assert_eq!(
        oriented
            .orientation()
            .oriented_cycle_from(IsogenyGraphNodeId(0)),
        Some(vec![
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(0)
        ])
    );
}
