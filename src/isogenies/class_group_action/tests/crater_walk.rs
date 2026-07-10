use super::{bu, cm_field_minus_23_curve, crater_report_with_nodes, split_three_ideal};
use crate::isogenies::{
    class_group_action::{CraterWalkReport, CraterWalkTermination},
    graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId, endomorphisms::CraterShape},
};

#[test]
fn crater_walk_report_records_a_closed_horizontal_cycle() {
    let graph = IsogenyGraphBuilder::new(cm_field_minus_23_curve(), 3)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_101 degree-three graph should build");

    let report = graph
        .crater_walk_report(split_three_ideal(), IsogenyGraphNodeId(0))
        .expect("walk report should build from crater evidence");

    assert_eq!(report.start(), IsogenyGraphNodeId(0));
    assert_eq!(report.ideal().norm(), &bu(3));
    assert_eq!(report.crater_shape(), CraterShape::Cycle { length: 3 });
    assert_eq!(
        report.visited(),
        &[
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2),
            IsogenyGraphNodeId(0)
        ]
    );
    assert_eq!(report.cycle_length(), Some(3));
    assert!(report.is_closed_cycle());
    assert!(report.start_in_crater());
    assert_eq!(report.termination(), CraterWalkTermination::ClosedCycle);
}

#[test]
fn crater_walk_report_records_non_crater_start_without_cycle() {
    let graph = IsogenyGraphBuilder::new(cm_field_minus_23_curve(), 3)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_101 degree-three graph should build");

    let report = graph
        .crater_walk_report(split_three_ideal(), IsogenyGraphNodeId(99))
        .expect("walk report should build even when the start is not in the crater");

    assert_eq!(report.start(), IsogenyGraphNodeId(99));
    assert_eq!(report.crater_shape(), CraterShape::Cycle { length: 3 });
    assert_eq!(report.visited(), &[IsogenyGraphNodeId(99)]);
    assert_eq!(report.cycle_length(), None);
    assert!(!report.is_closed_cycle());
    assert!(!report.start_in_crater());
    assert_eq!(
        report.termination(),
        CraterWalkTermination::StartOutsideCrater
    );
}

#[test]
fn crater_walk_report_distinguishes_crater_start_without_certified_outgoing_edge() {
    let crater = crater_report_with_nodes(bu(3), vec![IsogenyGraphNodeId(7)], Vec::new());

    let report =
        CraterWalkReport::from_crater_report(&crater, split_three_ideal(), IsogenyGraphNodeId(7));

    assert_eq!(report.start(), IsogenyGraphNodeId(7));
    assert_eq!(report.visited(), &[IsogenyGraphNodeId(7)]);
    assert_eq!(report.cycle_length(), None);
    assert!(!report.is_closed_cycle());
    assert!(report.start_in_crater());
    assert_eq!(
        report.termination(),
        CraterWalkTermination::NoCertifiedOutgoingEdge
    );
}
