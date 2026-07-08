use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::PrimeNormIdeal,
    quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
};
use crate::isogenies::{
    class_group_action::{HorizontalIdealReport, HorizontalIdealStatus},
    graphs::{
        IsogenyGraphEdgeId, IsogenyGraphNodeId,
        endomorphisms::{
            CraterReport, CraterShape, HorizontalEdgeReport, HorizontalEdgeStatus,
            VolcanoStructureReport,
        },
    },
};

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn order_minus_23() -> ImaginaryQuadraticOrder {
    ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(1))
        .expect("D = -23 should define an imaginary quadratic maximal order")
}

fn horizontal_edge(status: HorizontalEdgeStatus) -> HorizontalEdgeReport {
    HorizontalEdgeReport::new(
        IsogenyGraphEdgeId(7),
        IsogenyGraphNodeId(1),
        IsogenyGraphNodeId(2),
        status,
    )
}

fn crater_report(prime: BigUint, edges: Vec<HorizontalEdgeReport>) -> CraterReport {
    CraterReport::new(
        prime.clone(),
        VolcanoStructureReport::from_floor_paths(prime, Vec::new(), Vec::new()),
        Vec::new(),
        edges,
        CraterShape::EmptyCertifiedCrater,
    )
}

#[test]
fn certified_edge_with_matching_ideal_norm_yields_witness() {
    let edge = horizontal_edge(HorizontalEdgeStatus::CertifiedByAltitude);
    let ideal = PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23");

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
    let ideal = PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23");

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
    let ideal = PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23");

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
    let ideal = PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23");

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
    let ideal = PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23");

    let reports = HorizontalIdealReport::for_crater_report(&crater, ideal);

    assert_eq!(reports.len(), 1);
    assert_eq!(reports[0].status(), HorizontalIdealStatus::DegreeMismatch);
    assert!(reports[0].witness().is_none());
}
