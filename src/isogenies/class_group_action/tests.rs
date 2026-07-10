use num_bigint::BigUint;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::{
        binary_quadratic_forms::QuadraticClassGroup,
        quadratic_ideals::PrimeNormIdeal,
        quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
    },
};
use crate::isogenies::{
    class_group_action::{
        CraterIdealLabelError, CraterIdealLabelReport, CraterIdealPrimeBehavior,
        HorizontalIdealReport, HorizontalIdealStatus,
    },
    graphs::{
        IsogenyGraphBuilder, IsogenyGraphEdgeId, IsogenyGraphNodeId,
        endomorphisms::{
            CraterReport, CraterShape, HorizontalEdgeReport, HorizontalEdgeStatus,
            VolcanoStructureReport,
        },
    },
};

type F7 = crate::fields::Fp7;

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn order_minus_23() -> ImaginaryQuadraticOrder {
    ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(1))
        .expect("D = -23 should define an imaginary quadratic maximal order")
}

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid F_7 curve")
}

fn split_three_ideal() -> PrimeNormIdeal {
    PrimeNormIdeal::split(order_minus_23(), bu(3), bu(1))
        .expect("3 splits in the order of discriminant -23")
}

fn ramified_twenty_three_ideal() -> PrimeNormIdeal {
    PrimeNormIdeal::ramified(order_minus_23(), bu(23))
        .expect("23 ramifies in the order of discriminant -23")
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

fn class_group_minus_23() -> QuadraticClassGroup {
    QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should define an imaginary quadratic class group")
}

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
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");
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
        2
    );
}

#[test]
fn crater_ideal_label_accepts_split_ideal_for_matching_crater_and_class_group() {
    let crater = crater_report(bu(3), Vec::new());
    let ideal = split_three_ideal();

    let report = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect("split ideal should label the matching crater");

    assert_eq!(report.crater_prime(), &bu(3));
    assert_eq!(report.ideal().norm(), &bu(3));
    assert_eq!(report.prime_behavior(), CraterIdealPrimeBehavior::Split);
}

#[test]
fn crater_ideal_label_accepts_ramified_ideal_for_matching_crater_and_class_group() {
    let crater = crater_report(bu(23), Vec::new());
    let ideal = ramified_twenty_three_ideal();

    let report = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect("ramified ideal should label the matching crater");

    assert_eq!(report.crater_prime(), &bu(23));
    assert_eq!(report.ideal().norm(), &bu(23));
    assert_eq!(report.prime_behavior(), CraterIdealPrimeBehavior::Ramified);
}

#[test]
fn crater_ideal_label_rejects_crater_prime_mismatch() {
    let crater = crater_report(bu(5), Vec::new());
    let ideal = split_three_ideal();

    let error = CraterIdealLabelReport::new(&crater, &class_group_minus_23(), ideal)
        .expect_err("crater prime should match the ideal norm");

    assert_eq!(
        error,
        CraterIdealLabelError::PrimeNormMismatch {
            ideal_norm: bu(3),
            crater_prime: bu(5),
        }
    );
}

#[test]
fn crater_ideal_label_rejects_class_group_discriminant_mismatch() {
    let crater = crater_report(bu(3), Vec::new());
    let ideal = split_three_ideal();
    let wrong_class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-20))
        .expect("D = -20 should define an imaginary quadratic class group");

    let error = CraterIdealLabelReport::new(&crater, &wrong_class_group, ideal)
        .expect_err("ideal order and class group must have the same discriminant");

    assert_eq!(
        error,
        CraterIdealLabelError::OrderDiscriminantMismatch {
            ideal_discriminant: (-23).into(),
            class_group_discriminant: (-20).into(),
        }
    );
}

#[test]
fn crater_walk_report_records_a_closed_horizontal_cycle() {
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");

    let report = graph
        .crater_walk_report(split_three_ideal(), IsogenyGraphNodeId(0))
        .expect("walk report should build from crater evidence");

    assert_eq!(report.start(), IsogenyGraphNodeId(0));
    assert_eq!(report.ideal().norm(), &bu(3));
    assert_eq!(
        report.crater_shape(),
        CraterShape::TwoVertex {
            directed_edge_count: 2
        }
    );
    assert_eq!(
        report.visited(),
        &[
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(0)
        ]
    );
    assert_eq!(report.cycle_length(), Some(2));
    assert!(report.is_closed_cycle());
}

#[test]
fn crater_walk_report_records_non_crater_start_without_cycle() {
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");

    let report = graph
        .crater_walk_report(split_three_ideal(), IsogenyGraphNodeId(99))
        .expect("walk report should build even when the start is not in the crater");

    assert_eq!(report.start(), IsogenyGraphNodeId(99));
    assert_eq!(
        report.crater_shape(),
        CraterShape::TwoVertex {
            directed_edge_count: 2
        }
    );
    assert_eq!(report.visited(), &[IsogenyGraphNodeId(99)]);
    assert_eq!(report.cycle_length(), None);
    assert!(!report.is_closed_cycle());
}
