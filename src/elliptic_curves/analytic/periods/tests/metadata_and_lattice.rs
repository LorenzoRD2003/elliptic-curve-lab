use num_complex::Complex64;

use crate::elliptic_curves::analytic::lattice::HasAnalyticLatticeContext;
use crate::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve, ComplexLattice, CurvePeriodLatticeComparisonReport,
    LatticeSumTruncation, NumericalRecoveryMetadata, PeriodLatticeApprox, PeriodRecoveryMethod,
    PeriodRecoveryStatus, UpperHalfPlanePoint,
};
use crate::numerics::{ApproxTolerance, HasComplexApproxComparison};

#[test]
fn numerical_recovery_metadata_preserves_caller_supplied_fields() {
    let metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::AgmViaLegendre,
        PeriodRecoveryStatus::ValidationFailed,
        5,
        7,
        128,
        3,
        ApproxTolerance::strict(),
        Some(2.5e-10),
    );

    assert_eq!(
        metadata.resolved_method(),
        PeriodRecoveryMethod::AgmViaLegendre
    );
    assert_eq!(metadata.status(), PeriodRecoveryStatus::ValidationFailed);
    assert_eq!(metadata.newton_iterations_used(), 5);
    assert_eq!(metadata.agm_iterations_used(), 7);
    assert_eq!(metadata.integration_steps_used(), 128);
    assert_eq!(metadata.branch_lattice_searches_used(), 3);
    assert_eq!(metadata.tolerance(), ApproxTolerance::strict());
    assert_eq!(metadata.validation_residual_norm(), Some(2.5e-10));
    assert!(!metadata.succeeded());
}

#[test]
fn numerical_recovery_metadata_can_report_success_without_residual() {
    let metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::Hybrid,
        PeriodRecoveryStatus::Succeeded,
        4,
        6,
        64,
        1,
        ApproxTolerance::loose(),
        None,
    );

    assert_eq!(metadata.resolved_method(), PeriodRecoveryMethod::Hybrid);
    assert_eq!(metadata.status(), PeriodRecoveryStatus::Succeeded);
    assert_eq!(metadata.validation_residual_norm(), None);
    assert!(metadata.succeeded());
}

#[test]
fn standard_from_tau_uses_the_standard_z_plus_z_tau_basis() {
    let tau = UpperHalfPlanePoint::tau_i();
    let periods = PeriodLatticeApprox::standard_from_tau(tau.clone());

    assert_eq!(periods.omega1(), &Complex64::new(1.0, 0.0));
    assert_eq!(periods.omega2(), tau.tau());
    assert_eq!(periods.tau(), &tau);
}

#[test]
fn new_recovers_tau_from_the_supplied_lattice() {
    let lattice = ComplexLattice::new(Complex64::new(2.0, 0.0), Complex64::new(1.0, 2.0)).unwrap();
    let periods = PeriodLatticeApprox::new(lattice.clone()).unwrap();

    assert_eq!(periods.lattice(), &lattice);
    assert_eq!(periods.tau().tau(), &Complex64::new(0.5, 1.0));
}

#[test]
fn recovery_report_compares_recovered_and_curve_side_j_values() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let periods = PeriodLatticeApprox::standard_from_tau(tau);
    let recovered_j = curve.j_invariant().unwrap();

    let report = CurvePeriodLatticeComparisonReport::new(
        curve,
        periods,
        recovered_j,
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(report.recovered_j(), report.curve_j());
    assert_eq!(report.difference(), &Complex64::new(0.0, 0.0));
    assert!(report.agrees_approximately());
}

#[test]
fn recovery_report_reuses_the_shared_context_traits() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let periods = PeriodLatticeApprox::standard_from_tau(tau.clone());
    let report = CurvePeriodLatticeComparisonReport::new(
        curve,
        periods.clone(),
        Complex64::new(1728.0, 0.0),
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(report.tau(), periods.tau());
    assert_eq!(report.lattice(), periods.lattice());
    assert_eq!(report.left(), report.recovered_j());
    assert_eq!(report.right(), report.curve_j());
}
