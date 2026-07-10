use super::{bu, cm_field_minus_23_curve, crater_report, horizontal_edge, split_three_ideal};
use crate::isogenies::{
    class_group_action::{HorizontalIdealReport, HorizontalIdealStatus},
    graphs::{IsogenyGraphBuilder, IsogenyGraphEdgeId, endomorphisms::HorizontalEdgeStatus},
};

#[test]
fn certified_edge_with_matching_ideal_norm_yields_witness() {
    let edge = horizontal_edge(HorizontalEdgeStatus::CertifiedByAltitude);
    let ideal = split_three_ideal();

    let report = HorizontalIdealReport::from_certified_edge_and_ideal(edge, bu(3), ideal);

    assert_eq!(report.status(), HorizontalIdealStatus::CertifiedCompatible);
    let witness = report
        .witness()
        .expect("compatible report should carry a witness");
    assert_eq!(witness.edge().edge_id(), IsogenyGraphEdgeId(7));
    assert_eq!(witness.prime(), &bu(3));
    assert_eq!(witness.ideal().norm(), &bu(3));
}

#[test]
fn uncertified_edge_does_not_get_an_ideal_witness() {
    let edge = horizontal_edge(HorizontalEdgeStatus::SuspectedByWeakSurfaceEvidence);
    let ideal = split_three_ideal();

    let report = HorizontalIdealReport::from_certified_edge_and_ideal(edge, bu(3), ideal);

    assert_eq!(
        report.status(),
        HorizontalIdealStatus::EdgeNotCertifiedHorizontal
    );
    assert!(report.witness().is_none());
}

#[test]
fn degree_mismatch_is_reported_without_witness() {
    let edge = horizontal_edge(HorizontalEdgeStatus::CertifiedByAltitude);
    let ideal = split_three_ideal();

    let report = HorizontalIdealReport::from_certified_edge_and_ideal(edge, bu(5), ideal);

    assert_eq!(report.status(), HorizontalIdealStatus::DegreeMismatch);
    assert!(report.witness().is_none());
}

#[test]
fn crater_report_helper_uses_the_crater_prime_for_each_edge() {
    let crater = crater_report(
        bu(3),
        vec![
            horizontal_edge(HorizontalEdgeStatus::CertifiedByAltitude),
            horizontal_edge(HorizontalEdgeStatus::SuspectedByWeakSurfaceEvidence),
        ],
    );
    let ideal = split_three_ideal();

    let reports = HorizontalIdealReport::for_crater_report(&crater, ideal);

    assert_eq!(reports.len(), 2);
    assert_eq!(
        reports[0].status(),
        HorizontalIdealStatus::CertifiedCompatible
    );
    assert_eq!(
        reports[0]
            .witness()
            .expect("certified report should carry a witness")
            .prime(),
        &bu(3)
    );
    assert_eq!(
        reports[1].status(),
        HorizontalIdealStatus::EdgeNotCertifiedHorizontal
    );
}

#[test]
fn crater_report_helper_reports_mismatched_ideal_norms() {
    let crater = crater_report(
        bu(5),
        vec![horizontal_edge(HorizontalEdgeStatus::CertifiedByAltitude)],
    );
    let ideal = split_three_ideal();

    let reports = HorizontalIdealReport::for_crater_report(&crater, ideal);

    assert_eq!(reports.len(), 1);
    assert_eq!(reports[0].status(), HorizontalIdealStatus::DegreeMismatch);
    assert!(reports[0].witness().is_none());
}

#[test]
fn graph_method_uses_the_ideal_norm_as_the_crater_prime() {
    let graph = IsogenyGraphBuilder::new(cm_field_minus_23_curve(), 3)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_101 degree-three graph should build");
    let ideal = split_three_ideal();
    let crater = graph
        .volcano_crater_report(ideal.norm())
        .expect("crater report should build for the ideal norm");

    let direct_reports = HorizontalIdealReport::for_crater_report(&crater, ideal.clone());
    let method_reports = graph
        .horizontal_ideal_reports(ideal)
        .expect("graph method should reuse the ideal norm as the local volcano prime");

    assert_eq!(method_reports, direct_reports);
    assert_eq!(
        method_reports
            .iter()
            .filter(|report| report.status() == HorizontalIdealStatus::CertifiedCompatible)
            .count(),
        3
    );
}
