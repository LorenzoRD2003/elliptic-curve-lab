use num_complex::Complex64;
use proptest::prelude::*;

use crate::elliptic_curves::analytic::lattice::HasAnalyticLatticeContext;
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, CanonicalTauRecoveryReport,
    CompleteEllipticIntegralKApprox, CompleteEllipticIntegralKMetadata, ComplexAgmBranchChoice,
    ComplexAgmConfig, ComplexAgmStatus, ComplexLattice, CubicRootConfiguration,
    CubicRootConfigurationReport, CubicRootRecoveryReport, CubicRootSeparation,
    LatticeSumTruncation, LegendreOrbitElementKind, LegendreParameter,
    LegendreParameterConditioning, LegendreReduction, LegendreReductionReport,
    NumericalRecoveryMetadata, PeriodBasisRecoveryReport, PeriodLatticeApprox,
    PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryReport, PeriodRecoveryStatus,
    RecoveredPeriodBasis, RecoveredPeriodBasisReport, TauRecoveryReport, UpperHalfPlanePoint,
    WeierstrassCubicRoots, classify_cubic_root_configuration,
    classify_legendre_parameter_conditioning,
    complementary_complete_elliptic_integral_k_from_lambda,
    complete_elliptic_integral_k_from_lambda, complex_agm, complex_agm_trace,
    cubic_root_configuration_report, is_in_standard_fundamental_domain,
    legendre_period_integral_report, legendre_reduction_report, recover_canonical_tau_from_curve,
    recover_period_basis, recover_period_basis_from_legendre_reduction, recover_tau_from_curve,
    recover_weierstrass_cubic_roots, recover_weierstrass_cubic_roots_from_invariants,
    recover_weierstrass_cubic_roots_with_report,
};
use crate::fields::ComplexApprox;
use crate::numerics::HasComplexApproxComparison;

fn stable_real_split_curve_strategy() -> impl Strategy<Value = AnalyticWeierstrassCurve> {
    (0.4f64..3.0, 0.4f64..3.0)
        .prop_filter("real roots should stay well separated", |(e1, e2)| {
            let e3 = -(*e1 + *e2);
            (e1 - e2).abs() >= 0.2 && (e1 - e3).abs() >= 0.2 && (e2 - e3).abs() >= 0.2
        })
        .prop_map(|(e1, e2)| {
            let roots = WeierstrassCubicRoots::new(
                Complex64::new(e1, 0.0),
                Complex64::new(e2, 0.0),
                Complex64::new(-(e1 + e2), 0.0),
                ApproxTolerance::strict(),
            )
            .expect("strategy only yields distinct real roots");
            AnalyticWeierstrassCurve::new(roots.g2(), roots.g3())
                .expect("roots with distinct entries should define a nonsingular curve")
        })
}

#[test]
fn config_constructor_preserves_caller_supplied_values() {
    let tolerance = ApproxTolerance::new(1.0e-8, 2.0e-8);
    let config = PeriodRecoveryConfig::new(tolerance, 9, 7, 192, 3, 11).unwrap();

    assert_eq!(config.tolerance(), tolerance);
    assert_eq!(config.newton_max_iterations(), 9);
    assert_eq!(config.agm_max_iterations(), 7);
    assert_eq!(config.abel_jacobi_integration_steps(), 192);
    assert_eq!(config.branch_lattice_search_radius(), 3);
    assert_eq!(config.fundamental_domain_reduction_max_steps(), 11);
}

#[test]
fn config_rejects_zero_budgets() {
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 0, 1, 1, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 0, 1, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 0, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 1, 0, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 1, 1, 0),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
}

#[test]
fn config_presets_are_ordered_and_explicit() {
    let educational = PeriodRecoveryConfig::educational_default();
    let strict = PeriodRecoveryConfig::strict();
    let loose = PeriodRecoveryConfig::loose();

    assert_eq!(
        educational.tolerance(),
        ApproxTolerance::educational_default()
    );
    assert_eq!(educational.newton_max_iterations(), 12);
    assert_eq!(educational.agm_max_iterations(), 10);
    assert_eq!(educational.abel_jacobi_integration_steps(), 256);
    assert_eq!(educational.branch_lattice_search_radius(), 2);
    assert_eq!(educational.fundamental_domain_reduction_max_steps(), 16);

    assert_eq!(strict.tolerance(), ApproxTolerance::strict());
    assert!(strict.newton_max_iterations() > educational.newton_max_iterations());
    assert!(strict.agm_max_iterations() > educational.agm_max_iterations());
    assert!(strict.abel_jacobi_integration_steps() > educational.abel_jacobi_integration_steps());
    assert!(strict.branch_lattice_search_radius() > educational.branch_lattice_search_radius());
    assert!(
        strict.fundamental_domain_reduction_max_steps()
            > educational.fundamental_domain_reduction_max_steps()
    );

    assert_eq!(loose.tolerance(), ApproxTolerance::loose());
    assert!(loose.newton_max_iterations() < educational.newton_max_iterations());
    assert!(loose.agm_max_iterations() < educational.agm_max_iterations());
    assert!(loose.abel_jacobi_integration_steps() < educational.abel_jacobi_integration_steps());
    assert!(loose.branch_lattice_search_radius() < educational.branch_lattice_search_radius());
    assert!(
        loose.fundamental_domain_reduction_max_steps()
            < educational.fundamental_domain_reduction_max_steps()
    );
}

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
    let approximation =
        complete_elliptic_integral_k_from_lambda(&parameter, ComplexAgmConfig::strict()).unwrap();

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
    let approximation = complementary_complete_elliptic_integral_k_from_lambda(
        &parameter,
        ComplexAgmConfig::strict(),
    )
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
    let report =
        legendre_period_integral_report(&reduction, PeriodRecoveryConfig::strict()).unwrap();

    assert_eq!(report.lambda, reduction.parameter().clone());
    assert!(report.tau_candidate.re.is_finite());
    assert!(report.tau_candidate.im.is_finite());
    assert!(report.tau_candidate.im > 0.0);
    assert!(report.k_lambda.metadata().succeeded());
    assert!(report.k_complementary.metadata().succeeded());
}

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
    let integral_report =
        legendre_period_integral_report(&reduction, PeriodRecoveryConfig::strict()).unwrap();
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
    let report =
        recover_period_basis_from_legendre_reduction(&reduction, PeriodRecoveryConfig::strict())
            .unwrap();

    assert_eq!(report.reduction(), &reduction);
    assert!(ComplexApprox::eq_with_tolerance(
        report.basis().tau().tau(),
        &report.integral_report().tau_candidate,
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
    let k_report =
        legendre_period_integral_report(&reduction, PeriodRecoveryConfig::strict()).unwrap();
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
    let report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();

    assert_eq!(report.curve(), &curve);
    assert_eq!(report.legendre_reduction().roots(), report.roots());
    assert_eq!(report.tau(), report.periods().tau());
    assert!(ComplexApprox::eq_with_tolerance(
        report.tau().tau(),
        &report.k_report().tau_candidate,
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
    let k_report =
        legendre_period_integral_report(&reduction, PeriodRecoveryConfig::strict()).unwrap();
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
    let period_basis_report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();
    let tau_report = recover_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();

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
    let k_report =
        legendre_period_integral_report(&reduction, PeriodRecoveryConfig::strict()).unwrap();
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
        crate::elliptic_curves::reduce_tau_to_standard_fundamental_domain(tau_report.tau(), 8)
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
    let report = recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();

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
    let report = recover_canonical_tau_from_curve(&curve, one_step_config).unwrap();

    assert!(
        report.fundamental_domain_reduction().steps().len()
            <= one_step_config.fundamental_domain_reduction_max_steps()
    );
    assert_eq!(report.fundamental_domain_reduction().steps().len(), 1);
}

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
    assert_eq!(metadata.cardano_product_residual_norm(), None);
    assert_eq!(metadata.cardano_discriminant(), None);
    assert_eq!(metadata.selected_u_branch_index(), None);
    assert_eq!(metadata.selected_v_branch_index(), None);
    assert_eq!(metadata.used_principal_cardano_branches(), None);
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
fn numerical_recovery_metadata_can_attach_cardano_diagnostics() {
    let metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::Hybrid,
        PeriodRecoveryStatus::Succeeded,
        2,
        0,
        0,
        0,
        ApproxTolerance::strict(),
        Some(1.0e-12),
    )
    .with_cardano_diagnostics(Complex64::new(3.0, -4.0), 2.5e-14, 0, 2);

    assert_eq!(
        metadata.cardano_discriminant(),
        Some(&Complex64::new(3.0, -4.0))
    );
    assert_eq!(metadata.cardano_product_residual_norm(), Some(2.5e-14));
    assert_eq!(metadata.selected_u_branch_index(), Some(0));
    assert_eq!(metadata.selected_v_branch_index(), Some(2));
    assert_eq!(metadata.used_principal_cardano_branches(), Some(false));
}

#[test]
fn cubic_roots_preserve_input_order_without_claiming_a_canonical_sort() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(2.0, 0.5),
        Complex64::new(-3.0, 0.0),
        Complex64::new(1.0, -0.25),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let [first, second, third] = roots.roots();

    assert_eq!(first, &Complex64::new(2.0, 0.5));
    assert_eq!(second, &Complex64::new(-3.0, 0.0));
    assert_eq!(third, &Complex64::new(1.0, -0.25));
}

#[test]
fn cubic_roots_can_match_up_to_permutation() {
    let original = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let permuted = WeierstrassCubicRoots::new(
        Complex64::new(-3.0, 0.0),
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(
        original.matching_permutation(&permuted, ApproxTolerance::strict()),
        Some([1, 2, 0])
    );
    assert!(original.approximately_matches_up_to_permutation(&permuted, ApproxTolerance::strict()));
}

#[test]
fn cubic_roots_reject_approximately_repeated_entries() {
    assert_eq!(
        WeierstrassCubicRoots::new(
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0 + 5.0e-13, 0.0),
            Complex64::new(-2.0, 0.0),
            ApproxTolerance::strict(),
        ),
        Err(AnalyticCurveError::RepeatedCubicRoot)
    );
}

#[test]
fn cubic_roots_recover_the_expected_symmetric_sums_and_invariants() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(roots.sum(), Complex64::new(0.0, 0.0));
    assert_eq!(roots.pairwise_products_sum(), Complex64::new(-7.0, 0.0));
    assert_eq!(roots.product(), Complex64::new(-6.0, 0.0));
    assert_eq!(roots.x_squared_coefficient(), Complex64::new(0.0, 0.0));
    assert_eq!(roots.g2(), Complex64::new(28.0, 0.0));
    assert_eq!(roots.g3(), Complex64::new(-24.0, 0.0));
}

#[test]
fn cubic_roots_detect_depressed_and_non_depressed_cases() {
    let depressed = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let non_depressed = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(4.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert!(depressed.is_approximately_depressed(ApproxTolerance::strict()));
    assert!(!non_depressed.is_approximately_depressed(ApproxTolerance::strict()));
}

#[test]
fn cubic_roots_report_min_pairwise_distance() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(0.0, 0.0),
        Complex64::new(3.0, 4.0),
        Complex64::new(1.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(roots.min_pairwise_distance(), 1.0);
}

#[test]
fn legendre_parameter_new_rejects_exact_singular_and_non_finite_values() {
    assert_eq!(
        LegendreParameter::new(Complex64::new(0.0, 0.0)),
        Err(AnalyticCurveError::InvalidLegendreModulus)
    );
    assert_eq!(
        LegendreParameter::new(Complex64::new(1.0, 0.0)),
        Err(AnalyticCurveError::InvalidLegendreModulus)
    );
    assert_eq!(
        LegendreParameter::new(Complex64::new(f64::INFINITY, 0.0)),
        Err(AnalyticCurveError::InvalidLegendreModulus)
    );
}

#[test]
fn legendre_parameter_from_real_roots_chooses_a_deterministic_representative() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let parameter = LegendreParameter::from_roots(&roots, ApproxTolerance::strict()).unwrap();

    assert_eq!(parameter.lambda(), &Complex64::new(-0.25, 0.0));
    assert_eq!(parameter.one_minus_lambda(), Complex64::new(1.25, 0.0));
    assert!(!parameter.is_near_singular(ApproxTolerance::strict()));
}

#[test]
fn legendre_parameter_from_roots_is_permutation_invariant() {
    let original = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 1.0),
        Complex64::new(-2.0, 0.0),
        Complex64::new(1.0, -1.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let permuted = WeierstrassCubicRoots::new(
        Complex64::new(1.0, -1.0),
        Complex64::new(1.0, 1.0),
        Complex64::new(-2.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    let original_parameter =
        LegendreParameter::from_roots(&original, ApproxTolerance::strict()).unwrap();
    let permuted_parameter =
        LegendreParameter::from_roots(&permuted, ApproxTolerance::strict()).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        original_parameter.lambda(),
        permuted_parameter.lambda(),
        ApproxTolerance::strict()
    ));
}

#[test]
fn legendre_orbit_exposes_the_classical_six_transforms_in_fixed_order() {
    let parameter = LegendreParameter::new(Complex64::new(-0.25, 0.0)).unwrap();
    let orbit = parameter.orbit();
    let values = orbit.values();

    assert_eq!(
        orbit.element(LegendreOrbitElementKind::Lambda).lambda(),
        &Complex64::new(-0.25, 0.0)
    );
    assert_eq!(
        orbit
            .element(LegendreOrbitElementKind::OneMinusLambda)
            .lambda(),
        &Complex64::new(1.25, 0.0)
    );
    assert_eq!(
        orbit
            .element(LegendreOrbitElementKind::ReciprocalLambda)
            .lambda(),
        &Complex64::new(-4.0, 0.0)
    );
    assert_eq!(
        orbit
            .element(LegendreOrbitElementKind::ReciprocalOneMinusLambda)
            .lambda(),
        &Complex64::new(0.8, 0.0)
    );
    assert_eq!(
        orbit
            .element(LegendreOrbitElementKind::LambdaMinusOneOverLambda)
            .lambda(),
        &Complex64::new(5.0, 0.0)
    );
    assert_eq!(
        orbit
            .element(LegendreOrbitElementKind::LambdaOverLambdaMinusOne)
            .lambda(),
        &Complex64::new(0.2, 0.0)
    );
    assert_eq!(
        values,
        [
            Complex64::new(-0.25, 0.0),
            Complex64::new(1.25, 0.0),
            Complex64::new(-4.0, 0.0),
            Complex64::new(0.8, 0.0),
            Complex64::new(5.0, 0.0),
            Complex64::new(0.2, 0.0),
        ]
    );
}

#[test]
fn legendre_parameter_detects_near_zero_near_one_and_near_infinity() {
    let near_zero = LegendreParameter::new(Complex64::new(1.0e-13, 0.0)).unwrap();
    let near_one = LegendreParameter::new(Complex64::new(1.0 + 1.0e-13, 0.0)).unwrap();
    let near_infinity = LegendreParameter::new(Complex64::new(1.0e13, 0.0)).unwrap();
    let tolerance = ApproxTolerance::strict();

    assert!(near_zero.is_near_zero(tolerance));
    assert!(near_zero.is_near_singular(tolerance));
    assert!(near_one.is_near_one(tolerance));
    assert!(near_one.is_near_singular(tolerance));
    assert!(near_infinity.is_near_singular(tolerance));
}

#[test]
fn legendre_parameter_conditioning_classifies_all_three_singular_directions() {
    let tolerance = ApproxTolerance::strict();

    assert_eq!(
        classify_legendre_parameter_conditioning(
            &LegendreParameter::new(Complex64::new(-0.25, 0.0)).unwrap(),
            tolerance
        ),
        LegendreParameterConditioning::Generic
    );
    assert_eq!(
        classify_legendre_parameter_conditioning(
            &LegendreParameter::new(Complex64::new(1.0e-13, 0.0)).unwrap(),
            tolerance
        ),
        LegendreParameterConditioning::NearZero
    );
    assert_eq!(
        classify_legendre_parameter_conditioning(
            &LegendreParameter::new(Complex64::new(1.0 + 1.0e-13, 0.0)).unwrap(),
            tolerance
        ),
        LegendreParameterConditioning::NearOne
    );
    assert_eq!(
        classify_legendre_parameter_conditioning(
            &LegendreParameter::new(Complex64::new(1.0e13, 0.0)).unwrap(),
            tolerance
        ),
        LegendreParameterConditioning::NearInfinity
    );
    assert!(!LegendreParameterConditioning::Generic.is_near_singular());
    assert!(LegendreParameterConditioning::NearInfinity.is_near_singular());
}

#[test]
fn legendre_reduction_records_the_selected_affine_normalization() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let selected = reduction.selected_root_triple();

    assert_eq!(reduction.parameter().lambda(), &Complex64::new(-0.25, 0.0));
    assert_eq!(reduction.selected_permutation(), [2, 0, 1]);
    assert_eq!(selected[0], &Complex64::new(-3.0, 0.0));
    assert_eq!(selected[1], &Complex64::new(1.0, 0.0));
    assert_eq!(selected[2], &Complex64::new(2.0, 0.0));
    assert_eq!(reduction.x_translation(), Complex64::new(1.0, 0.0));
    assert_eq!(reduction.x_scale(), Complex64::new(-4.0, 0.0));
    assert_eq!(
        reduction.legendre_rhs_scale_factor(),
        Complex64::new(-256.0, 0.0)
    );
    assert_eq!(reduction.principal_sqrt_x_scale(), Complex64::new(0.0, 2.0));
    assert_eq!(reduction.legendre_y_scale(), Complex64::new(0.0, -16.0));
    assert_eq!(
        reduction.invariant_differential_scale(),
        Complex64::new(0.0, -0.25)
    );
}

#[test]
fn legendre_reduction_maps_chosen_roots_to_zero_one_and_lambda() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let selected = reduction.selected_root_triple();

    assert_eq!(
        reduction.legendre_x_from_original_x(*selected[0]),
        Complex64::new(1.0, 0.0)
    );
    assert_eq!(
        reduction.legendre_x_from_original_x(*selected[1]),
        Complex64::new(0.0, 0.0)
    );
    assert_eq!(
        reduction.legendre_x_from_original_x(*selected[2]),
        *reduction.parameter().lambda()
    );
    assert_eq!(
        reduction.original_x_from_legendre_x(Complex64::new(0.0, 0.0)),
        Complex64::new(1.0, 0.0)
    );
    assert_eq!(
        reduction.original_x_from_legendre_x(Complex64::new(1.0, 0.0)),
        Complex64::new(-3.0, 0.0)
    );
}

#[test]
fn legendre_reduction_reconstructs_the_original_cubic_from_legendre_x() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let legendre_x = Complex64::new(0.5, 0.0);
    let original_x = reduction.original_x_from_legendre_x(legendre_x);
    let direct_original_cubic = Complex64::new(4.0, 0.0)
        * (original_x - Complex64::new(1.0, 0.0))
        * (original_x - Complex64::new(2.0, 0.0))
        * (original_x - Complex64::new(-3.0, 0.0));

    assert_eq!(
        reduction.evaluate_original_cubic_from_legendre_x(legendre_x),
        direct_original_cubic
    );
}

#[test]
fn legendre_reduction_principal_branch_scales_are_algebraically_consistent() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
    let alpha = reduction.principal_sqrt_x_scale();
    let y_scale = reduction.legendre_y_scale();
    let differential_scale = reduction.invariant_differential_scale();

    assert_eq!(alpha.powu(2), reduction.x_scale());
    assert_eq!(y_scale.powu(2), reduction.legendre_rhs_scale_factor());
    assert_eq!(differential_scale * y_scale, reduction.x_scale());
}

#[test]
fn legendre_reduction_report_wraps_the_reduction_and_records_conditioning() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let report = LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap();

    assert_eq!(report.parameter().lambda(), &Complex64::new(-0.25, 0.0));
    assert_eq!(
        report.selected_orbit_element_relative_to_input_order(),
        LegendreOrbitElementKind::ReciprocalOneMinusLambda
    );
    assert_eq!(
        report.conditioning(),
        LegendreParameterConditioning::Generic
    );
    assert_eq!(report.tolerance(), ApproxTolerance::strict());
    assert!(ApproxTolerance::strict().real_close(report.singularity_distance(), 0.25));
    assert!(!report.is_near_singular());
    assert_eq!(report.reduction().selected_permutation(), [2, 0, 1]);
}

#[test]
fn legendre_reduction_report_helper_delegates_to_report_constructor() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(
        legendre_reduction_report(&roots, ApproxTolerance::strict()).unwrap(),
        LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap()
    );
}

#[test]
fn cubic_root_configuration_detects_three_approximately_real_roots() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(-3.0, 1.0e-13),
        Complex64::new(1.0, -2.0e-13),
        Complex64::new(2.0, 1.5e-13),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let report = cubic_root_configuration_report(&roots, ApproxTolerance::strict());

    assert_eq!(
        classify_cubic_root_configuration(&roots, ApproxTolerance::strict()),
        CubicRootConfiguration::ThreeApproximatelyReal
    );
    assert_eq!(
        report.configuration(),
        CubicRootConfiguration::ThreeApproximatelyReal
    );
    assert_eq!(report.separation(), CubicRootSeparation::WellSeparated);
    assert_eq!(report.conjugate_pair_residual(), None);
}

#[test]
fn cubic_root_configuration_detects_one_real_plus_conjugate_pair() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(2.0, 1.0),
        Complex64::new(-3.0, 0.0),
        Complex64::new(2.0, -1.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let report = CubicRootConfigurationReport::new(roots.clone(), ApproxTolerance::strict());

    assert_eq!(
        report.configuration(),
        CubicRootConfiguration::OneApproximatelyRealTwoApproximatelyConjugate
    );
    assert_eq!(report.separation(), CubicRootSeparation::WellSeparated);
    assert_eq!(report.tolerance(), ApproxTolerance::strict());
    assert!(report.conjugate_pair_residual().is_some());
    assert!(report.conjugate_pair_residual().unwrap() <= ApproxTolerance::strict().absolute);
    assert_eq!(report.roots(), &roots);
}

#[test]
fn cubic_root_configuration_detects_generic_complex_case() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 1.0),
        Complex64::new(-0.2, 0.1),
        Complex64::new(-0.8, -1.1),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let report = cubic_root_configuration_report(&roots, ApproxTolerance::strict());

    assert_eq!(
        report.configuration(),
        CubicRootConfiguration::GenericComplex
    );
    assert_eq!(report.separation(), CubicRootSeparation::WellSeparated);
    assert_eq!(report.conjugate_pair_residual(), None);
}

#[test]
fn cubic_root_configuration_report_tracks_nearly_repeated_status_under_looser_tolerance() {
    let roots = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(1.0 + 5.0e-8, 0.0),
        Complex64::new(-2.0 - 5.0e-8, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let report = cubic_root_configuration_report(&roots, ApproxTolerance::loose());

    assert_eq!(
        report.configuration(),
        CubicRootConfiguration::ThreeApproximatelyReal
    );
    assert_eq!(report.separation(), CubicRootSeparation::NearlyRepeated);
    assert!(ApproxTolerance::strict().real_close(report.min_pairwise_distance(), 5.0e-8));
}

#[test]
fn recovered_cubic_roots_match_a_real_split_example() {
    let expected = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(28.0, 0.0), Complex64::new(-24.0, 0.0))
            .unwrap();
    let roots = recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()).unwrap();

    assert!(roots.is_approximately_depressed(ApproxTolerance::strict()));
    assert!(roots.approximately_matches_up_to_permutation(&expected, ApproxTolerance::strict()));
    assert!(ComplexApprox::eq_with_tolerance(
        &roots.g2(),
        curve.g2(),
        ApproxTolerance::strict()
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        &roots.g3(),
        curve.g3(),
        ApproxTolerance::strict()
    ));
}

#[test]
fn recovered_cubic_roots_match_a_generic_complex_example_via_symmetric_invariants() {
    let source = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 1.0),
        Complex64::new(-2.0, 0.0),
        Complex64::new(1.0, -1.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(source.g2(), source.g3()).unwrap();
    let recovered =
        recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()).unwrap();

    assert!(recovered.is_approximately_depressed(ApproxTolerance::strict()));
    assert!(recovered.approximately_matches_up_to_permutation(&source, ApproxTolerance::strict()));
    assert!(ComplexApprox::eq_with_tolerance(
        &recovered.g2(),
        &source.g2(),
        ApproxTolerance::strict()
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        &recovered.g3(),
        &source.g3(),
        ApproxTolerance::strict()
    ));
}

#[test]
fn curve_and_invariants_recovery_surfaces_agree() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(28.0, 0.0), Complex64::new(-24.0, 0.0))
            .unwrap();
    let from_curve =
        recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()).unwrap();
    let from_invariants = recover_weierstrass_cubic_roots_from_invariants(
        curve.g2(),
        curve.g3(),
        PeriodRecoveryConfig::strict(),
    )
    .unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        &from_curve.g2(),
        &from_invariants.g2(),
        ApproxTolerance::strict()
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        &from_curve.g3(),
        &from_invariants.g3(),
        ApproxTolerance::strict()
    ));
}

#[test]
fn cubic_root_recovery_report_reuses_coefficient_comparisons_and_metadata() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(28.0, 0.0), Complex64::new(-24.0, 0.0))
            .unwrap();
    let report: CubicRootRecoveryReport =
        recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
            .unwrap();

    assert_eq!(report.curve(), &curve);
    assert!(report.g2_comparison().agrees_approximately());
    assert!(report.g3_comparison().agrees_approximately());
    assert!(report.reconstruction_agrees());
    assert_eq!(report.curve_g2(), curve.g2());
    assert_eq!(report.curve_g3(), curve.g3());
    assert_eq!(report.metadata().status(), PeriodRecoveryStatus::Succeeded);
    assert_eq!(
        report.metadata().resolved_method(),
        PeriodRecoveryMethod::Hybrid
    );
    assert_eq!(report.metadata().agm_iterations_used(), 0);
    assert_eq!(report.metadata().integration_steps_used(), 0);
    assert_eq!(report.metadata().branch_lattice_searches_used(), 0);
    assert!(report.metadata().validation_residual_norm().is_some());
    assert!(report.metadata().cardano_discriminant().is_some());
    assert!(report.metadata().cardano_product_residual_norm().is_some());
    assert!(report.metadata().selected_u_branch_index().is_some());
    assert!(report.metadata().selected_v_branch_index().is_some());
}

#[test]
fn ill_conditioned_complex_example_uses_newton_polishing() {
    let source = WeierstrassCubicRoots::new(
        Complex64::new(8.813789020059971, -6.296193572032816),
        Complex64::new(-5.70258988712044, -4.026550473696494),
        Complex64::new(-3.1111991329395314, 10.32274404572931),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(source.g2(), source.g3()).unwrap();
    let report =
        recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
            .unwrap();

    assert!(report.metadata().newton_iterations_used() > 0);
    assert!(report.metadata().cardano_product_residual_norm().is_some());
    assert!(report.metadata().cardano_discriminant().is_some());
    assert!(report.metadata().selected_u_branch_index().is_some());
    assert!(report.metadata().selected_v_branch_index().is_some());
    assert!(report.reconstruction_agrees());
    assert!(
        report
            .roots()
            .approximately_matches_up_to_permutation(&source, ApproxTolerance::strict())
    );
}

#[test]
fn slightly_perturbed_invariants_still_recover_a_consistent_real_configuration() {
    let source = WeierstrassCubicRoots::new(
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
        ApproxTolerance::strict(),
    )
    .unwrap();
    let perturbed_curve = AnalyticWeierstrassCurve::new(
        source.g2() + Complex64::new(1.0e-8, -3.0e-9),
        source.g3() + Complex64::new(-2.0e-8, 4.0e-9),
    )
    .unwrap();
    let report = recover_weierstrass_cubic_roots_with_report(
        &perturbed_curve,
        PeriodRecoveryConfig::strict(),
    )
    .unwrap();
    let strict_classification =
        cubic_root_configuration_report(report.roots(), report.metadata().tolerance());
    let loose_classification =
        cubic_root_configuration_report(report.roots(), ApproxTolerance::loose());

    assert!(report.reconstruction_agrees());
    assert!(
        report
            .roots()
            .approximately_matches_up_to_permutation(&source, ApproxTolerance::loose())
    );
    assert_eq!(
        strict_classification.configuration(),
        CubicRootConfiguration::GenericComplex
    );
    assert_eq!(
        loose_classification.configuration(),
        CubicRootConfiguration::ThreeApproximatelyReal
    );
    assert_eq!(
        loose_classification.separation(),
        CubicRootSeparation::WellSeparated
    );
}

#[test]
fn very_large_complex_invariants_can_trigger_branch_choice_ambiguity() {
    let source = WeierstrassCubicRoots::new(
        Complex64::new(6697.015551473084, 104627.82280425371),
        Complex64::new(91139.91685017172, -55531.26160040997),
        Complex64::new(-97836.93240164481, -49096.56120384374),
        ApproxTolerance::loose(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(source.g2(), source.g3()).unwrap();

    assert_eq!(
        recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()),
        Err(AnalyticCurveError::BranchChoiceAmbiguous)
    );
}

#[test]
fn very_large_complex_invariants_can_trigger_recovery_failure() {
    let source = WeierstrassCubicRoots::new(
        Complex64::new(-539492.1494784288, 538219.9587749569),
        Complex64::new(411534.7845088351, 160515.8843419563),
        Complex64::new(127957.36496959365, -698735.8431169132),
        ApproxTolerance::loose(),
    )
    .unwrap();
    let curve = AnalyticWeierstrassCurve::new(source.g2(), source.g3()).unwrap();

    assert_eq!(
        recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()),
        Err(AnalyticCurveError::CubicRootRecoveryFailed)
    );
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

    let report =
        PeriodRecoveryReport::new(curve, periods, recovered_j, ApproxTolerance::strict()).unwrap();

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
    let report = PeriodRecoveryReport::new(
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

#[test]
fn recovering_periods_from_a_curve_built_from_tau_recovers_the_full_lattice_scale() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();
    let source = ComplexLattice::from_tau(tau);

    assert!(
        ApproxTolerance::new(1.0e-2, 1.0e-2)
            .real_close(report.periods().covolume(), source.covolume())
    );
}

#[test]
fn recovering_tau_rho_from_a_hexagonal_lattice_curve_recovers_the_expected_modular_class() {
    let tau = UpperHalfPlanePoint::tau_rho();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report = recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        report.canonical_tau().tau(),
        tau.tau(),
        ApproxTolerance::new(1.0e-3, 1.0e-3)
    ));
}

#[test]
fn near_pure_cubic_root_recovery_handles_tau_rho_with_strict_tolerance() {
    let tau = UpperHalfPlanePoint::tau_rho();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report =
        recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
            .unwrap();

    assert!(report.reconstruction_agrees());
    assert!(report.roots().min_pairwise_distance() > 1.0e-3);
    assert_eq!(report.metadata().cardano_discriminant(), None);
    assert_eq!(report.metadata().cardano_product_residual_norm(), None);
    assert_eq!(report.metadata().selected_u_branch_index(), None);
    assert_eq!(report.metadata().selected_v_branch_index(), None);
}

#[test]
fn all_root_orderings_produce_modularly_equivalent_recovered_taus() {
    let roots = [
        Complex64::new(1.0, 0.0),
        Complex64::new(2.0, 0.0),
        Complex64::new(-3.0, 0.0),
    ];
    let permutations = [
        [0usize, 1usize, 2usize],
        [0usize, 2usize, 1usize],
        [1usize, 0usize, 2usize],
        [1usize, 2usize, 0usize],
        [2usize, 0usize, 1usize],
        [2usize, 1usize, 0usize],
    ];

    let canonical_taus = permutations
        .into_iter()
        .map(|indices| {
            let ordered = WeierstrassCubicRoots::new(
                roots[indices[0]],
                roots[indices[1]],
                roots[indices[2]],
                ApproxTolerance::strict(),
            )
            .unwrap();
            let curve = AnalyticWeierstrassCurve::new(ordered.g2(), ordered.g3()).unwrap();
            recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict())
                .unwrap()
                .canonical_tau()
                .clone()
        })
        .collect::<Vec<_>>();

    for tau in canonical_taus.iter().skip(1) {
        assert!(ComplexApprox::eq_with_tolerance(
            tau.tau(),
            canonical_taus[0].tau(),
            ApproxTolerance::new(1.0e-3, 1.0e-3)
        ));
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn canonical_tau_recovery_for_stable_real_split_curves_produces_a_valid_fundamental_domain_representative(
        curve in stable_real_split_curve_strategy(),
    ) {
        let config = PeriodRecoveryConfig::strict();
        let report = recover_canonical_tau_from_curve(&curve, config)
            .expect("stable real-split test family should recover a canonical tau");

        prop_assert!(report.metadata().succeeded());
        prop_assert!(report.periods().covolume().is_sign_positive());
        prop_assert!(is_in_standard_fundamental_domain(
            report.canonical_tau(),
            config.tolerance(),
        ));
        prop_assert_eq!(report.original_tau(), report.tau_recovery_report().tau());
        prop_assert_eq!(
            report.canonical_tau(),
            report.fundamental_domain_reduction().reduced_tau()
        );

        let matrix_image = report
            .accumulated_matrix()
            .apply(&report.original_tau())
            .expect("accumulated modular matrix should act on the recovered tau");

        prop_assert!(ComplexApprox::eq_with_tolerance(
            matrix_image.tau(),
            report.canonical_tau().tau(),
            config.tolerance(),
        ));
        prop_assert!(
            report.fundamental_domain_reduction().steps().len()
                <= config.fundamental_domain_reduction_max_steps()
        );
        prop_assert!(report.fundamental_domain_reduction().is_reduced());
    }
}
