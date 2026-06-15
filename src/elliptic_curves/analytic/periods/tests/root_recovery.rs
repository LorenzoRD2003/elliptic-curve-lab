use super::*;

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
    let roots = curve
        .recover_weierstrass_cubic_roots(PeriodRecoveryConfig::strict())
        .unwrap();

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
    let recovered = curve
        .recover_weierstrass_cubic_roots(PeriodRecoveryConfig::strict())
        .unwrap();

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
    let from_curve = curve
        .recover_weierstrass_cubic_roots(PeriodRecoveryConfig::strict())
        .unwrap();
    let from_invariants = WeierstrassCubicRoots::from_invariants(
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
    let report: CubicRootRecoveryReport = curve
        .recover_weierstrass_cubic_roots_with_report(PeriodRecoveryConfig::strict())
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
    assert!(report.cardano_discriminant().is_some());
    assert!(report.cardano_product_residual_norm().is_some());
    assert!(report.selected_u_branch_index().is_some());
    assert!(report.selected_v_branch_index().is_some());
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
    let report = curve
        .recover_weierstrass_cubic_roots_with_report(PeriodRecoveryConfig::strict())
        .unwrap();

    assert!(report.metadata().newton_iterations_used() > 0);
    assert!(report.cardano_product_residual_norm().is_some());
    assert!(report.cardano_discriminant().is_some());
    assert!(report.selected_u_branch_index().is_some());
    assert!(report.selected_v_branch_index().is_some());
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
    let report = perturbed_curve
        .recover_weierstrass_cubic_roots_with_report(PeriodRecoveryConfig::strict())
        .unwrap();
    let strict_classification = report
        .roots()
        .configuration_report(report.metadata().tolerance());
    let loose_classification = report
        .roots()
        .configuration_report(ApproxTolerance::loose());

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
        curve.recover_weierstrass_cubic_roots(PeriodRecoveryConfig::strict()),
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
        curve.recover_weierstrass_cubic_roots(PeriodRecoveryConfig::strict()),
        Err(AnalyticCurveError::CubicRootRecoveryFailed)
    );
}

#[test]
fn near_pure_cubic_root_recovery_handles_tau_rho_with_strict_tolerance() {
    let tau = UpperHalfPlanePoint::tau_rho();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report = curve
        .recover_weierstrass_cubic_roots_with_report(PeriodRecoveryConfig::strict())
        .unwrap();

    assert!(report.reconstruction_agrees());
    assert!(report.roots().min_pairwise_distance() > 1.0e-3);
    assert_eq!(report.cardano_discriminant(), None);
    assert_eq!(report.cardano_product_residual_norm(), None);
    assert_eq!(report.selected_u_branch_index(), None);
    assert_eq!(report.selected_v_branch_index(), None);
}
