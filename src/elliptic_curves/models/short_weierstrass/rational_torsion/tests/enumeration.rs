use num_bigint::BigInt;

use crate::elliptic_curves::{
    AffinePoint,
    short_weierstrass::rational_torsion::{
        enumeration::LutzNagellCandidateReport, integral_model::RationalIntegralModel,
    },
    traits::CurveModel,
};

use super::fixtures::{
    cyclic_six_fixture, product_two_two_fixture, q, rational_scaled_fixture,
    trivial_torsion_fixture,
};

#[test]
fn lutz_nagell_candidates_for_full_two_torsion_are_the_three_affine_roots_plus_identity() {
    let fixture = product_two_two_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("fixture should already be integral");
    let report = LutzNagellCandidateReport::from_integral_model(&witness)
        .expect("candidate enumeration should succeed");

    assert_eq!(report.discriminant(), &BigInt::from(64));
    assert_eq!(report.candidate_count(), 4);
    assert_eq!(
        report.points(),
        &[
            AffinePoint::infinity(),
            AffinePoint::new(q(-1, 1), q(0, 1)),
            AffinePoint::new(q(0, 1), q(0, 1)),
            AffinePoint::new(q(1, 1), q(0, 1)),
        ]
    );
}

#[test]
fn lutz_nagell_candidates_include_odd_order_pairs() {
    let fixture = cyclic_six_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("fixture should already be integral");
    let report = LutzNagellCandidateReport::from_integral_model(&witness)
        .expect("candidate enumeration should succeed");

    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(0, 1), q(1, 1)))
    );
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(0, 1), q(-1, 1)))
    );
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
    for point in report.points() {
        assert!(witness.curve().contains(point));
    }
}

#[test]
fn lutz_nagell_candidates_run_on_the_scaled_integral_companion() {
    let fixture = rational_scaled_fixture();
    let witness =
        RationalIntegralModel::from_curve(fixture.curve).expect("fixture should scale integrally");
    let report = LutzNagellCandidateReport::from_integral_model(&witness)
        .expect("candidate enumeration should succeed");

    assert_eq!(witness.scale(), &q(2, 1));
    assert_eq!(report.discriminant(), &BigInt::from(64));
    assert_eq!(report.candidate_count(), 4);
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(-1, 1), q(0, 1)))
    );
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(0, 1), q(0, 1)))
    );
    assert!(
        report
            .points()
            .contains(&AffinePoint::new(q(1, 1), q(0, 1)))
    );
}

#[test]
fn lutz_nagell_candidates_validate_every_returned_point_exactly() {
    let fixture = trivial_torsion_fixture();
    let witness = RationalIntegralModel::from_curve(fixture.curve)
        .expect("fixture should already be integral");
    let report = LutzNagellCandidateReport::from_integral_model(&witness)
        .expect("candidate enumeration should succeed");

    assert!(report.points().contains(&AffinePoint::infinity()));
    for point in report.points() {
        assert!(witness.curve().contains(point));
    }
}
