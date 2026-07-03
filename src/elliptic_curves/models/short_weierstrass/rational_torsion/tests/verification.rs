use crate::elliptic_curves::{
    AffinePoint,
    short_weierstrass::rational_torsion::{
        RationalTorsionGroupShape, RationalTorsionReport, integral_model::RationalIntegralModel,
    },
    traits::{CurveModel, GroupCurveModel},
};

use super::fixtures::{
    cyclic_five_fixture, cyclic_seven_fixture, cyclic_six_fixture, product_two_two_fixture, q,
    rational_scaled_fixture, trivial_torsion_fixture,
};

#[test]
fn rational_torsion_report_classifies_full_two_torsion() {
    let fixture = product_two_two_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve.clone())
        .expect("fixture should already be integral");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }
    );
    assert_eq!(report.points(), fixture.sample_points.as_slice());
    assert_eq!(report.candidate_count(), 4);
    assert_eq!(report.rejected_candidate_count(), 0);
}

#[test]
fn rational_torsion_report_classifies_cyclic_six() {
    let fixture = cyclic_six_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve.clone())
        .expect("fixture should already be integral");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::Cyclic { order: 6 }
    );
    assert_eq!(report.points().len(), 6);
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(2, 1), q(3, 1)))
    );
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(2, 1), q(-3, 1)))
    );
}

#[test]
fn rational_torsion_report_classifies_cyclic_five_and_seven() {
    for (fixture, order, generator) in [
        (
            cyclic_five_fixture(),
            5,
            AffinePoint::new(q(-1, 3), q(1, 2)),
        ),
        (
            cyclic_seven_fixture(),
            7,
            AffinePoint::new(q(3, 1), q(8, 1)),
        ),
    ] {
        let report = fixture
            .curve
            .rational_torsion()
            .expect("cyclic torsion fixture should classify");

        assert_eq!(
            report.group().shape(),
            RationalTorsionGroupShape::Cyclic { order }
        );
        assert_eq!(report.points().len(), order);
        assert!(report.points().contains(&generator));
        assert_eq!(
            fixture
                .curve
                .exact_mazur_order(&generator)
                .expect("generator should be on the curve"),
            Some(order)
        );
    }
}

#[test]
fn rational_torsion_report_classifies_trivial_torsion() {
    let fixture = trivial_torsion_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("fixture should already be integral");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    assert_eq!(report.group().shape(), RationalTorsionGroupShape::Trivial);
    assert_eq!(report.points(), &[AffinePoint::infinity()]);
    assert!(report.candidate_count() >= report.points().len());
}

#[test]
fn lutz_nagell_candidates_can_be_rejected_as_nontorsion() {
    let fixture = trivial_torsion_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("fixture should already be integral");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    assert_eq!(report.group().shape(), RationalTorsionGroupShape::Trivial);
    assert_eq!(report.points(), &[AffinePoint::infinity()]);
    assert_eq!(report.candidate_count(), 3);
    assert_eq!(report.rejected_candidate_count(), 2);
}

#[test]
fn rational_torsion_report_transports_scaled_points_back_to_source_curve() {
    let fixture = rational_scaled_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve.clone())
        .expect("fixture should scale integrally");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    assert_eq!(
        report.group().shape(),
        RationalTorsionGroupShape::ProductZ2Z2m { m: 1 }
    );
    assert_eq!(report.points(), fixture.sample_points.as_slice());
    for point in report.points() {
        assert!(fixture.curve.contains(point));
    }
}

#[test]
fn reported_point_orders_are_exact() {
    let fixture = cyclic_six_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve.clone())
        .expect("fixture should already be integral");
    let report =
        RationalTorsionReport::from_integral_model(&witness).expect("torsion should classify");

    for point in report.points() {
        let order = fixture
            .curve
            .exact_mazur_order(point)
            .expect("reported point should have a well-defined order")
            .expect("reported point should be torsion");
        assert!(
            fixture.curve.is_identity(
                &fixture
                    .curve
                    .mul_scalar(point, order)
                    .expect("scalar multiple should stay on curve")
            )
        );
        for proper_divisor in 1..order {
            if order.is_multiple_of(proper_divisor) {
                assert!(
                    !fixture.curve.is_identity(
                        &fixture
                            .curve
                            .mul_scalar(point, proper_divisor)
                            .expect("scalar multiple should stay on curve")
                    )
                );
            }
        }
    }
}
