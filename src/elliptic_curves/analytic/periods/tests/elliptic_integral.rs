use super::*;
use num_complex::Complex64;

#[test]
fn complete_elliptic_integral_k_approx_preserves_caller_supplied_fields() {
    let parameter = LegendreParameter::new(Complex64::new(0.25, -0.5)).unwrap();
    let metadata = CompleteEllipticIntegralKMetadata::new(
        PeriodRecoveryMethod::AgmViaLegendre,
        PeriodRecoveryStatus::Succeeded,
        ApproxTolerance::strict(),
        7,
        Complex64::new(0.75, 0.0).sqrt(),
        true,
    );
    let approximation = CompleteEllipticIntegralKApprox::new(
        parameter.clone(),
        Complex64::new(1.854074677, 0.125),
        metadata.clone(),
    );

    assert_eq!(approximation.parameter(), &parameter);
    assert_eq!(approximation.value(), &Complex64::new(1.854074677, 0.125));
    assert_eq!(approximation.metadata(), &metadata);
    assert_eq!(
        approximation.metadata().resolved_method(),
        PeriodRecoveryMethod::AgmViaLegendre
    );
    assert_eq!(
        approximation.metadata().status(),
        PeriodRecoveryStatus::Succeeded
    );
    assert_eq!(
        approximation.metadata().tolerance(),
        ApproxTolerance::strict()
    );
    assert_eq!(approximation.metadata().agm_iterations_used(), 7);
    assert!(
        approximation
            .metadata()
            .used_principal_complementary_branch()
    );
    assert!(approximation.metadata().succeeded());
}

#[test]
fn complete_elliptic_integral_k_from_lambda_matches_the_known_half_parameter_value() {
    let parameter = LegendreParameter::new(Complex64::new(0.5, 0.0)).unwrap();
    let approximation = parameter
        .complete_elliptic_integral_k(ComplexAgmConfig::strict())
        .unwrap();

    assert_eq!(approximation.parameter(), &parameter);
    assert_eq!(
        approximation.metadata().resolved_method(),
        PeriodRecoveryMethod::AgmViaLegendre
    );
    assert!(ComplexApprox::eq_with_tolerance(
        approximation.value(),
        &Complex64::new(1.854_074_677_301_371_9, 0.0),
        ApproxTolerance::new(1.0e-12, 1.0e-12)
    ));
}

#[test]
fn complementary_complete_elliptic_integral_from_half_parameter_matches_the_same_value() {
    let parameter = LegendreParameter::new(Complex64::new(0.5, 0.0)).unwrap();
    let approximation = parameter
        .complementary_complete_elliptic_integral_k(ComplexAgmConfig::strict())
        .unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        approximation.value(),
        &Complex64::new(1.854_074_677_301_371_9, 0.0),
        ApproxTolerance::new(1.0e-12, 1.0e-12)
    ));
}

#[test]
fn legendre_period_integral_report_recovers_tau_i_for_lambda_one_half() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.5, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let report = reduction
        .period_integral_report(PeriodRecoveryConfig::strict())
        .unwrap();

    assert_eq!(report.lambda(), reduction.parameter());
    assert!(report.tau_candidate().re.is_finite());
    assert!(report.tau_candidate().im.is_finite());
    assert!(report.tau_candidate().im > 0.0);
    assert!(report.k_lambda().metadata().succeeded());
    assert!(report.k_complementary().metadata().succeeded());
}
