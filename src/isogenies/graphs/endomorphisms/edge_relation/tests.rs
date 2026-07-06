use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::QuadraticDiscriminant,
};
use crate::isogenies::graphs::endomorphisms::{
    IsogenyEdgeEndomorphismReport, IsogenyEdgeEndomorphismTentativeRelation,
    edge_relation::infer_relation,
};

fn candidate_set(discriminant: i64) -> EndomorphismRingCandidateSet {
    QuadraticDiscriminant::new(discriminant)
        .factorization()
        .expect("test discriminant should factor canonically")
        .endomorphism_ring_candidates()
        .expect("candidate orders should construct")
}

#[test]
fn identical_singleton_levels_are_possibly_horizontal() {
    let source = candidate_set(-4);
    let target = candidate_set(-4);

    let report =
        IsogenyEdgeEndomorphismReport::from_candidate_sets(&BigUint::from(2u8), &source, &target)
            .expect("report should build");

    assert_eq!(report.source_possible_levels(), &[0]);
    assert_eq!(report.target_possible_levels(), &[0]);
    assert_eq!(
        report.relation(),
        &IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
    );
}

#[test]
fn different_quadratic_fields_are_unsupported() {
    let source = candidate_set(-16);
    let target = candidate_set(-3);

    let report = source
        .tentative_edge_endomorphism_report(&BigUint::from(2u8), &target)
        .expect("report should build");

    assert_eq!(report.source_possible_levels(), &[0, 1]);
    assert_eq!(report.target_possible_levels(), &[0]);
    assert_eq!(
        report.relation(),
        &IsogenyEdgeEndomorphismTentativeRelation::Unsupported
    );
}

#[test]
fn mixed_horizontal_and_vertical_possibilities_are_ambiguous() {
    let source = candidate_set(-16);
    let target = candidate_set(-16);

    let report =
        IsogenyEdgeEndomorphismReport::from_candidate_sets(&BigUint::from(2u8), &source, &target)
            .expect("report should build");

    assert_eq!(report.source_possible_levels(), &[0, 1]);
    assert_eq!(report.target_possible_levels(), &[0, 1]);
    assert_eq!(
        report.relation(),
        &IsogenyEdgeEndomorphismTentativeRelation::Ambiguous
    );
}

#[test]
fn level_classifier_recovers_all_tentative_variants() {
    assert_eq!(
        infer_relation(&[0], &[0]),
        IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal
    );
    assert_eq!(
        infer_relation(&[1], &[0]),
        IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending
    );
    assert_eq!(
        infer_relation(&[0], &[1]),
        IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending
    );
    assert_eq!(
        infer_relation(&[0, 1], &[0]),
        IsogenyEdgeEndomorphismTentativeRelation::Ambiguous
    );
    assert_eq!(
        infer_relation(&[0], &[2]),
        IsogenyEdgeEndomorphismTentativeRelation::Unsupported
    );
}

#[test]
fn relation_can_be_classified_from_floor_distances() {
    assert_eq!(
        IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(2, 2),
        Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyHorizontal)
    );
    assert_eq!(
        IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(2, 1),
        Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyDescending)
    );
    assert_eq!(
        IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(1, 2),
        Some(IsogenyEdgeEndomorphismTentativeRelation::PossiblyAscending)
    );
    assert_eq!(
        IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(3, 1),
        None
    );
}
