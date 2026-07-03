use crate::elliptic_curves::short_weierstrass::rational_torsion::{
    RationalTorsionGroupShape, RationalTorsionReport,
};
use crate::elliptic_curves::traits::CurveModel;

use super::fixtures::{
    cyclic_six_fixture, product_two_two_fixture, q, rational_scaled_fixture,
    trivial_torsion_fixture,
};

#[test]
fn stage_zero_fixtures_are_valid_curves_with_documented_sample_points() {
    for fixture in [
        product_two_two_fixture(),
        cyclic_six_fixture(),
        trivial_torsion_fixture(),
        rational_scaled_fixture(),
    ] {
        for point in &fixture.sample_points {
            assert!(
                fixture.curve.contains(point),
                "{} should contain its documented sample point {point:?}",
                fixture.name
            );
        }
    }
}

#[test]
fn stage_zero_report_keeps_points_as_the_canonical_payload() {
    let fixture = product_two_two_fixture();
    let report = RationalTorsionReport::new(
        fixture.curve.clone(),
        fixture.curve,
        q(1, 1),
        fixture.expected_group,
        fixture.sample_points.clone(),
        fixture.sample_points.len(),
        0,
    );

    assert_eq!(report.original_curve(), report.integral_model());
    assert_eq!(report.scale(), &q(1, 1));
    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }
    );
    assert_eq!(report.points(), fixture.sample_points.as_slice());
    assert_eq!(report.candidate_count(), fixture.sample_points.len());
    assert_eq!(report.rejected_candidate_count(), 0);
}
