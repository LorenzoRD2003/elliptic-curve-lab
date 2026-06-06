use crate::elliptic_curves::analytic::{
    ComplexAnalyticCurveLabReport, EllipticFunctionTruncation, LatticeSumTruncation, SpecialJKind,
    SpecialTauKind, TorusToCurveValues, UniformizationExperimentReport, UpperHalfPlanePoint,
};
use crate::fields::ComplexApprox;
use num_complex::Complex64;

#[test]
fn lab_report_keeps_tau_lattice_and_q_parameter_consistent() {
    let tau = UpperHalfPlanePoint::tau_generic_example();
    let report = ComplexAnalyticCurveLabReport::from_tau(
        tau.clone(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert_eq!(report.tau(), &tau);
    assert_eq!(
        report.lattice(),
        &crate::elliptic_curves::analytic::ComplexLattice::from_tau(tau.clone())
    );
    assert_eq!(report.q_parameter().tau(), &tau);
    assert_eq!(report.lattice().tau().unwrap(), tau);
}

#[test]
fn lab_report_keeps_j_consistent_across_all_model_surfaces() {
    let report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_generic_example(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        &report.invariants().j_invariant,
        &report.analytic_curve().j_invariant().unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        &report.invariants().j_invariant,
        report.short_model().j_invariant(),
        crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
    ));
}

#[test]
fn special_tau_kind_distinguishes_i_rho_and_generic_cases() {
    let tau_i_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_i(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let tau_rho_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_rho(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let generic_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_generic_example(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert_eq!(tau_i_report.special_tau_kind(), SpecialTauKind::TauI);
    assert_eq!(tau_rho_report.special_tau_kind(), SpecialTauKind::TauRho);
    assert_eq!(generic_report.special_tau_kind(), SpecialTauKind::Generic);
}

#[test]
fn special_j_kind_reflects_the_classical_special_values() {
    let tau_i_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_i(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let tau_rho_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_rho(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let generic_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_generic_example(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert_eq!(tau_i_report.special_j_kind(), SpecialJKind::Near1728);
    assert_eq!(tau_rho_report.special_j_kind(), SpecialJKind::NearZero);
    assert_eq!(generic_report.special_j_kind(), SpecialJKind::Generic);
}

#[test]
fn uniformization_report_derives_global_curve_membership_from_samples() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.3, 0.2),
            Complex64::new(0.5, 0.0),
        ],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert_eq!(
        report.all_points_lie_on_curve(),
        report
            .sampled_points()
            .iter()
            .all(|point| point.lies_on_curve())
    );
}

#[test]
fn uniformization_report_keeps_the_same_curve_across_all_samples() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![Complex64::new(0.3, 0.2), Complex64::new(0.5, 0.0)],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert!(
        report
            .sampled_points()
            .iter()
            .all(|point| point.curve() == report.curve())
    );
}

#[test]
fn uniformization_report_can_include_both_finite_points_and_a_pole() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![Complex64::new(0.0, 0.0), Complex64::new(0.3, 0.2)],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert!(matches!(
        report.sampled_points()[0].values(),
        TorusToCurveValues::Pole
    ));
    assert!(matches!(
        report.sampled_points()[1].values(),
        TorusToCurveValues::FiniteValues { .. }
    ));
}
