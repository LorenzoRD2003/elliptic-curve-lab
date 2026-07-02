use super::*;
use num_complex::Complex64;

#[test]
fn recovered_period_basis_wraps_one_validated_lattice() {
    let basis =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(1.0, 3.0)).unwrap();

    assert_eq!(basis.omega1(), &Complex64::new(2.0, 0.0));
    assert_eq!(basis.omega2(), &Complex64::new(1.0, 3.0));
    assert_eq!(basis.tau().tau(), &Complex64::new(0.5, 1.5));
    assert_eq!(basis.oriented_area(), 6.0);
    assert_eq!(basis.covolume(), 6.0);
}

#[test]
fn recovered_period_basis_report_preserves_caller_supplied_fields() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let integral_report = reduction
        .period_integral_report(PeriodRecoveryConfig::strict())
        .unwrap();
    let basis =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();
    let report =
        RecoveredPeriodBasisReport::new(reduction.clone(), integral_report.clone(), basis.clone());

    assert_eq!(report.reduction(), &reduction);
    assert_eq!(report.integral_report(), &integral_report);
    assert_eq!(
        report.invariant_differential_scale(),
        reduction.invariant_differential_scale()
    );
    assert_eq!(report.basis(), &basis);
}

#[test]
fn recover_period_basis_from_legendre_reduction_builds_a_valid_basis() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let report = reduction
        .recover_period_basis(PeriodRecoveryConfig::strict())
        .unwrap();

    assert_eq!(report.reduction(), &reduction);
    assert!(ComplexApprox::eq_with_tolerance(
        report.basis().tau().tau(),
        report.integral_report().tau_candidate(),
        ApproxTolerance::strict()
    ));
    assert!(report.basis().oriented_area().is_sign_positive());
}

#[test]
fn period_basis_recovery_report_preserves_caller_supplied_fields() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)).unwrap();
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let k_report = reduction
        .period_integral_report(PeriodRecoveryConfig::strict())
        .unwrap();
    let periods =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();
    let tau = periods.tau();
    let basis_report =
        RecoveredPeriodBasisReport::new(reduction.clone(), k_report.clone(), periods.clone());
    let metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::Hybrid,
        PeriodRecoveryStatus::Succeeded,
        4,
        6,
        0,
        0,
        ApproxTolerance::strict(),
        Some(1.0e-12),
    );
    let report = PeriodBasisRecoveryReport::new(
        curve.clone(),
        roots.clone(),
        basis_report.clone(),
        metadata.clone(),
    );

    assert_eq!(report.curve(), &curve);
    assert_eq!(report.roots(), &roots);
    assert_eq!(report.basis_report(), &basis_report);
    assert_eq!(report.legendre_reduction(), &reduction);
    assert_eq!(report.k_report(), &k_report);
    assert_eq!(report.periods(), &periods);
    assert_eq!(report.tau(), tau);
    assert_eq!(report.metadata(), &metadata);
}

#[test]
fn recover_period_basis_builds_a_complete_curve_level_report() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
    let report = curve
        .recover_period_basis(PeriodRecoveryConfig::strict())
        .unwrap();

    assert_eq!(report.curve(), &curve);
    assert_eq!(report.legendre_reduction().roots(), report.roots());
    assert_eq!(report.tau(), report.periods().tau());
    assert!(ComplexApprox::eq_with_tolerance(
        report.tau().tau(),
        report.k_report().tau_candidate(),
        ApproxTolerance::strict()
    ));
    assert_eq!(
        report.metadata().resolved_method(),
        PeriodRecoveryMethod::Hybrid
    );
    assert!(
        report.metadata().newton_iterations_used()
            <= 3 * PeriodRecoveryConfig::strict().newton_max_iterations()
    );
    assert!(report.metadata().agm_iterations_used() > 0);
    assert!(report.metadata().succeeded());
}

#[test]
fn tau_recovery_report_preserves_caller_supplied_fields() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)).unwrap();
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let k_report = reduction
        .period_integral_report(PeriodRecoveryConfig::strict())
        .unwrap();
    let periods =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();
    let basis_report =
        RecoveredPeriodBasisReport::new(reduction.clone(), k_report.clone(), periods.clone());
    let period_basis_report = PeriodBasisRecoveryReport::new(
        curve.clone(),
        roots.clone(),
        basis_report.clone(),
        NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::Succeeded,
            4,
            6,
            0,
            0,
            ApproxTolerance::strict(),
            Some(1.0e-12),
        ),
    );
    let report = TauRecoveryReport::new(period_basis_report.clone());

    assert_eq!(report.period_basis_report(), &period_basis_report);
    assert_eq!(report.curve(), &curve);
    assert_eq!(report.roots(), &roots);
    assert_eq!(report.basis_report(), &basis_report);
    assert_eq!(report.legendre_reduction(), &reduction);
    assert_eq!(report.k_report(), &k_report);
    assert_eq!(report.periods(), &periods);
    assert_eq!(report.tau(), periods.tau());
    assert_eq!(report.metadata(), period_basis_report.metadata());
}

#[test]
fn recover_tau_from_curve_reuses_period_basis_recovery() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
    let period_basis_report = curve
        .recover_period_basis(PeriodRecoveryConfig::strict())
        .unwrap();
    let tau_report = curve.recover_tau(PeriodRecoveryConfig::strict()).unwrap();

    assert_eq!(tau_report.period_basis_report(), &period_basis_report);
    assert_eq!(tau_report.tau(), period_basis_report.tau());
    assert_eq!(tau_report.periods(), period_basis_report.periods());
    assert_eq!(tau_report.metadata(), period_basis_report.metadata());
}

#[test]
fn canonical_tau_recovery_report_preserves_caller_supplied_fields() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)).unwrap();
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let k_report = reduction
        .period_integral_report(PeriodRecoveryConfig::strict())
        .unwrap();
    let periods =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();
    let basis_report =
        RecoveredPeriodBasisReport::new(reduction.clone(), k_report.clone(), periods.clone());
    let period_basis_report = PeriodBasisRecoveryReport::new(
        curve.clone(),
        roots.clone(),
        basis_report.clone(),
        NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::Succeeded,
            4,
            6,
            0,
            0,
            ApproxTolerance::strict(),
            Some(1.0e-12),
        ),
    );
    let tau_report = TauRecoveryReport::new(period_basis_report.clone());
    let fundamental_domain_reduction =
        crate::elliptic_curves::analytic::reduce_tau_to_standard_fundamental_domain(
            tau_report.tau(),
            8,
        )
        .unwrap();
    let report =
        CanonicalTauRecoveryReport::new(tau_report.clone(), fundamental_domain_reduction.clone());

    assert_eq!(report.tau_recovery_report(), &tau_report);
    assert_eq!(
        report.fundamental_domain_reduction(),
        &fundamental_domain_reduction
    );
    assert_eq!(report.curve(), &curve);
    assert_eq!(report.roots(), &roots);
    assert_eq!(report.basis_report(), &basis_report);
    assert_eq!(report.periods(), &periods);
    assert_eq!(report.original_tau(), tau_report.tau());
    assert_eq!(
        report.canonical_tau(),
        fundamental_domain_reduction.reduced_tau()
    );
    assert_eq!(
        report.accumulated_matrix(),
        fundamental_domain_reduction.accumulated_matrix()
    );
    assert_eq!(report.metadata(), period_basis_report.metadata());
}

#[test]
fn recover_canonical_tau_from_curve_reduces_the_natural_tau_to_the_standard_domain() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
    let report = curve
        .recover_canonical_tau(PeriodRecoveryConfig::strict())
        .unwrap();

    assert!(is_in_standard_fundamental_domain(
        report.canonical_tau(),
        ApproxTolerance::strict()
    ));
    assert_eq!(report.original_tau(), report.tau_recovery_report().tau());
    assert_eq!(
        report.canonical_tau(),
        report.fundamental_domain_reduction().reduced_tau()
    );
    let matrix_image = report
        .accumulated_matrix()
        .apply(&report.original_tau())
        .unwrap();
    assert!(ComplexApprox::eq_with_tolerance(
        matrix_image.tau(),
        report.canonical_tau().tau(),
        ApproxTolerance::strict()
    ));
}

#[test]
fn recover_canonical_tau_from_curve_uses_the_configured_modular_reduction_budget() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
    let one_step_config =
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 20, 16, 512, 4, 1).unwrap();
    let report = curve.recover_canonical_tau(one_step_config).unwrap();

    assert!(
        report.fundamental_domain_reduction().steps().len()
            <= one_step_config.fundamental_domain_reduction_max_steps()
    );
    assert_eq!(report.fundamental_domain_reduction().steps().len(), 1);
}
