use super::*;

#[test]
fn complex_agm_config_constructor_preserves_caller_supplied_values() {
    let tolerance = ApproxTolerance::new(1.0e-8, 2.0e-8);
    let config = ComplexAgmConfig::new(tolerance, 9).unwrap();

    assert_eq!(config.tolerance(), tolerance);
    assert_eq!(config.max_iterations(), 9);
}

#[test]
fn complex_agm_config_rejects_invalid_tolerances_and_zero_budget() {
    assert_eq!(
        ComplexAgmConfig::new(ApproxTolerance::new(0.0, 0.0), 4),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        ComplexAgmConfig::new(ApproxTolerance::new(-1.0, 1.0e-9), 4),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        ComplexAgmConfig::new(ApproxTolerance::new(1.0e-9, f64::NAN), 4),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        ComplexAgmConfig::new(ApproxTolerance::strict(), 0),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
}

#[test]
fn complex_agm_config_can_be_derived_from_period_recovery_config() {
    let config = ComplexAgmConfig::from_period_recovery_config(PeriodRecoveryConfig::strict());

    assert_eq!(config.tolerance(), ApproxTolerance::strict());
    assert_eq!(
        config.max_iterations(),
        PeriodRecoveryConfig::strict().agm_max_iterations()
    );
}

#[test]
fn complex_agm_presets_match_the_period_recovery_agm_presets() {
    assert_eq!(
        ComplexAgmConfig::educational_default(),
        ComplexAgmConfig::from_period_recovery_config(PeriodRecoveryConfig::educational_default())
    );
    assert_eq!(
        ComplexAgmConfig::strict(),
        ComplexAgmConfig::from_period_recovery_config(PeriodRecoveryConfig::strict())
    );
    assert_eq!(
        ComplexAgmConfig::loose(),
        ComplexAgmConfig::from_period_recovery_config(PeriodRecoveryConfig::loose())
    );
}

#[test]
fn complex_agm_succeeds_immediately_for_equal_inputs() {
    let value = Complex64::new(1.25, -0.5);
    let result = complex_agm(value, value, ComplexAgmConfig::strict()).unwrap();
    let trace = complex_agm_trace(value, value, ComplexAgmConfig::strict()).unwrap();

    assert_eq!(result.status(), ComplexAgmStatus::Succeeded);
    assert_eq!(result.iterations_used(), 0);
    assert_eq!(result.input_a(), &value);
    assert_eq!(result.input_b(), &value);
    assert_eq!(result.final_a(), &value);
    assert_eq!(result.final_b(), &value);
    assert_eq!(result.agm(), &value);
    assert_eq!(result.final_gap_norm(), 0.0);
    assert!(result.succeeded());

    assert_eq!(trace.iterations(), &[]);
    assert_eq!(trace.result(), &result);
}

#[test]
fn complex_agm_trace_records_the_selected_square_root_branch() {
    let config = ComplexAgmConfig::new(ApproxTolerance::new(1.0, 1.0), 1).unwrap();
    let trace =
        complex_agm_trace(Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0), config).unwrap();
    let [step] = trace.iterations() else {
        panic!("expected exactly one recorded AGM step");
    };

    let expected_next_a = Complex64::new(0.5, 0.5);
    let expected_principal = Complex64::new(
        std::f64::consts::FRAC_1_SQRT_2,
        std::f64::consts::FRAC_1_SQRT_2,
    );

    assert_eq!(step.index(), 0);
    assert_eq!(step.a_n(), &Complex64::new(1.0, 0.0));
    assert_eq!(step.b_n(), &Complex64::new(0.0, 1.0));
    assert!(ComplexApprox::eq_with_tolerance(
        step.next_a(),
        &expected_next_a,
        ApproxTolerance::strict()
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        step.principal_sqrt_product(),
        &expected_principal,
        ApproxTolerance::strict()
    ));
    assert_eq!(
        step.selected_branch(),
        ComplexAgmBranchChoice::PrincipalSqrt
    );
    assert_eq!(step.selected_geometric_mean(), step.next_b());
    assert!(step.next_gap_norm() < 0.5);
}

#[test]
fn complex_agm_can_hit_the_iteration_limit_without_erroring() {
    let config = ComplexAgmConfig::new(ApproxTolerance::strict(), 1).unwrap();
    let result = complex_agm(Complex64::new(1.0, 0.0), Complex64::new(0.5, 0.0), config).unwrap();

    assert_eq!(result.status(), ComplexAgmStatus::HitIterationLimit);
    assert_eq!(result.iterations_used(), 1);
    assert!(result.final_gap_norm().is_finite());
    assert!(result.final_gap_norm() < 0.5);
    assert!(!result.succeeded());
}

#[test]
fn complex_agm_positive_real_inputs_converge_to_a_common_limit() {
    let result = complex_agm(
        Complex64::new(1.0, 0.0),
        Complex64::new(0.5, 0.0),
        ComplexAgmConfig::strict(),
    )
    .unwrap();

    assert_eq!(result.status(), ComplexAgmStatus::Succeeded);
    assert!(result.iterations_used() > 0);
    assert!(ComplexApprox::eq_with_tolerance(
        result.final_a(),
        result.final_b(),
        ComplexAgmConfig::strict().tolerance()
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        result.agm(),
        &Complex64::new(0.7283955155234534, 0.0),
        ApproxTolerance::new(1.0e-12, 1.0e-12)
    ));
}

#[test]
fn complex_agm_rejects_non_finite_inputs() {
    assert_eq!(
        complex_agm(
            Complex64::new(f64::INFINITY, 0.0),
            Complex64::new(1.0, 0.0),
            ComplexAgmConfig::strict(),
        ),
        Err(AnalyticCurveError::InvalidAgmInput)
    );
    assert_eq!(
        complex_agm_trace(
            Complex64::new(1.0, 0.0),
            Complex64::new(f64::NAN, 0.0),
            ComplexAgmConfig::strict(),
        ),
        Err(AnalyticCurveError::InvalidAgmInput)
    );
}
