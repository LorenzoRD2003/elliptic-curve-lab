use num_complex::Complex64;

use super::{
    AbelJacobiConfig, AbelJacobiContourReport, AbelJacobiInitialBranchChoice,
    AbelJacobiIntegralApprox, AbelJacobiIntegralDecomposition, AbelJacobiIntegralNumerics,
    AbelJacobiPointRecoveryReport, AbelJacobiRecoveryMetadata, AbelJacobiRecoveryStatus,
    AbelJacobiRoundtripValidationReport, AbelJacobiValidationPolicy,
    InvariantRecoveryInterpretation, InverseUniformizationPointRecoveryReport,
    LegendreContourStrategy, PointRoundTripValidationConfig, PointRoundTripValidationReport,
    approximate_abel_jacobi_integral, recover_torus_point_from_curve_point,
    recover_torus_point_from_curve_point_with_periods,
    validate_canonical_tau_recovery_by_j_invariant,
    validate_point_inverse_uniformization_roundtrip,
    validate_point_inverse_uniformization_roundtrip_with_periods,
    validate_recovered_lattice_invariants, validate_recovered_tau_by_j_invariant,
    validate_tau_recovery_report_by_j_invariant,
};
use crate::elliptic_curves::analytic::periods::{
    PeriodRecoveryConfig, RecoveredPeriodBasis, recover_canonical_tau_from_curve,
    recover_period_basis, recover_tau_from_curve, recover_weierstrass_cubic_roots,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, AnalyticWeierstrassCurve, ApproxTolerance,
    ComplexApproxComparison, ComplexLattice, EllipticFunctionTruncation, HasAnalyticLatticeContext,
    HasComplexApproxComparison, LatticeSumTruncation, UpperHalfPlanePoint,
    map_torus_point_to_curve,
};

#[test]
fn inverse_uniformization_validation_report_recomputes_tau_side_invariants_and_j() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();

    let report =
        validate_recovered_tau_by_j_invariant(&curve, &tau, truncation, ApproxTolerance::strict())
            .unwrap();

    assert_eq!(report.curve(), &curve);
    assert_eq!(report.tau(), &tau);
    assert_eq!(report.lattice(), &ComplexLattice::from_tau(tau.clone()));
    assert_eq!(report.lattice_truncation(), truncation);
    assert_eq!(
        report.recovered_j(),
        &report.recovered_invariants().j_invariant
    );
    assert_eq!(report.curve_j(), &curve.j_invariant().unwrap());
    assert_eq!(report.difference(), &Complex64::new(0.0, 0.0));
    assert!(report.agrees_approximately());
}

#[test]
fn inverse_uniformization_validation_report_reuses_shared_context_traits() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report =
        validate_recovered_tau_by_j_invariant(&curve, &tau, truncation, ApproxTolerance::loose())
            .unwrap();

    assert_eq!(report.tau(), &tau);
    assert_eq!(report.lattice(), &ComplexLattice::from_tau(tau.clone()));
    assert_eq!(report.left(), report.recovered_j());
    assert_eq!(report.right(), report.curve_j());
}

#[test]
fn tau_recovery_validation_wrapper_uses_the_natural_recovered_tau() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let tau_report = recover_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();

    let validation = validate_tau_recovery_report_by_j_invariant(
        &tau_report,
        truncation,
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(validation.curve(), tau_report.curve());
    assert_eq!(validation.tau(), &tau_report.tau());
    assert!(validation.recovered_j().re.is_finite());
    assert!(validation.recovered_j().im.is_finite());
}

#[test]
fn canonical_tau_validation_wrapper_uses_the_canonical_tau() {
    let tau = UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let canonical_report =
        recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();

    let validation = validate_canonical_tau_recovery_by_j_invariant(
        &canonical_report,
        truncation,
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(validation.curve(), canonical_report.curve());
    assert_eq!(validation.tau(), canonical_report.canonical_tau());
    assert!(validation.recovered_j().re.is_finite());
    assert!(validation.recovered_j().im.is_finite());
}

#[test]
fn invariant_recovery_validation_detects_direct_agreement_for_matching_lattice() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let periods = RecoveredPeriodBasis::from_lattice(ComplexLattice::from_tau(tau.clone()));

    let report = validate_recovered_lattice_invariants(
        &curve,
        &periods,
        truncation,
        ApproxTolerance::strict(),
    )
    .unwrap();

    assert_eq!(
        report.interpretation(),
        InvariantRecoveryInterpretation::DirectAgreement
    );
    assert!(report.direct_scale_sensitive_agreement());
    assert!(report.same_j_invariant_approximately());
    assert_eq!(report.tau(), &tau);
    assert_eq!(report.lattice(), periods.lattice());
}

#[test]
fn invariant_recovery_validation_detects_same_j_but_scale_sensitive_mismatch() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let periods =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();

    let report = validate_recovered_lattice_invariants(
        &curve,
        &periods,
        truncation,
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(
        report.interpretation(),
        InvariantRecoveryInterpretation::SameModularClassButScaleSensitiveMismatch
    );
    assert!(!report.direct_scale_sensitive_agreement());
    assert!(report.same_j_invariant_approximately());
}

#[test]
fn invariant_recovery_validation_detects_inconsistent_recovery() {
    let tau_curve = UpperHalfPlanePoint::tau_i();
    let tau_wrong = UpperHalfPlanePoint::from_re_im(0.3, 1.2).unwrap();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau_curve, truncation).unwrap();
    let periods = RecoveredPeriodBasis::from_lattice(ComplexLattice::from_tau(tau_wrong));

    let report = validate_recovered_lattice_invariants(
        &curve,
        &periods,
        truncation,
        ApproxTolerance::loose(),
    )
    .unwrap();

    assert_eq!(
        report.interpretation(),
        InvariantRecoveryInterpretation::Inconsistent
    );
    assert!(!report.same_j_invariant_approximately());
}

#[test]
fn abel_jacobi_config_presets_stay_positive() {
    for config in [
        AbelJacobiConfig::educational_default(),
        AbelJacobiConfig::strict(),
        AbelJacobiConfig::loose(),
    ] {
        assert_eq!(
            config.legendre_contour_strategy,
            LegendreContourStrategy::CanonicalSegmentThenRay
        );
        assert!(config.integration_steps > 0);
        assert!(config.segment_samples > 0);
        assert!(config.ray_samples > 0);
        assert!(config.max_branch_adjustments > 0);
        assert!(config.max_lattice_corrections > 0);
        assert!(config.validation_policy.lattice_truncation_radius > 0);
        assert!(config.validation_policy.function_truncation_radius > 0);
    }
}

#[test]
fn point_roundtrip_validation_config_presets_stay_positive() {
    for config in [
        PointRoundTripValidationConfig::educational_default(),
        PointRoundTripValidationConfig::strict(),
        PointRoundTripValidationConfig::loose(),
    ] {
        assert!(config.lattice_truncation().radius() > 0);
        assert!(config.function_truncation().radius() > 0);
    }
}

#[test]
fn abel_jacobi_contour_report_preserves_geometric_choices() {
    let contour = AbelJacobiContourReport::new(
        LegendreContourStrategy::CanonicalSegmentThenRay,
        Complex64::new(0.5, -0.25),
        Complex64::new(3.0, 4.0),
        1.2,
        5.0,
        20.0,
        0.8,
    )
    .unwrap();

    assert_eq!(
        contour.legendre_contour_strategy(),
        LegendreContourStrategy::CanonicalSegmentThenRay
    );
    assert_eq!(contour.start(), &Complex64::new(0.5, -0.25));
    assert_eq!(contour.anchor(), &Complex64::new(3.0, 4.0));
    assert_eq!(contour.theta(), 1.2);
    assert_eq!(contour.radius(), 5.0);
    assert_eq!(contour.tail_length(), 20.0);
    assert_eq!(contour.min_distance_to_branch_points(), 0.8);
}

#[test]
fn abel_jacobi_integral_approx_preserves_explicit_fields() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let point = AnalyticCurvePoint::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0));
    let contour = AbelJacobiContourReport::new(
        LegendreContourStrategy::CanonicalSegmentThenRay,
        Complex64::new(1.0, 2.0),
        Complex64::new(3.0, 4.0),
        0.25,
        5.0,
        20.0,
        0.75,
    )
    .unwrap();
    let approx = AbelJacobiIntegralApprox::new(
        curve.clone(),
        point.clone(),
        contour.clone(),
        Complex64::new(0.25, -0.5),
        AbelJacobiIntegralDecomposition::new(
            AbelJacobiInitialBranchChoice::Alternate,
            Complex64::new(0.1, 0.2),
            Complex64::new(0.05, -0.1),
            Complex64::new(0.1, -0.6),
        )
        .unwrap(),
        AbelJacobiIntegralNumerics::new(128, 3, ApproxTolerance::strict()),
    )
    .unwrap();

    assert_eq!(approx.curve(), &curve);
    assert_eq!(approx.point(), &point);
    assert_eq!(approx.contour(), &contour);
    assert_eq!(approx.value(), &Complex64::new(0.25, -0.5));
    assert_eq!(
        approx.initial_branch_choice(),
        AbelJacobiInitialBranchChoice::Alternate
    );
    assert_eq!(approx.segment_integral(), &Complex64::new(0.1, 0.2));
    assert_eq!(approx.ray_integral(), &Complex64::new(0.05, -0.1));
    assert_eq!(approx.tail_correction(), &Complex64::new(0.1, -0.6));
    assert_eq!(approx.integration_steps_used(), 128);
    assert_eq!(approx.branch_adjustments_used(), 3);
    assert_eq!(approx.tolerance(), ApproxTolerance::strict());
}

#[test]
fn approximate_abel_jacobi_integral_at_infinity_returns_exact_zero() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let point = AnalyticCurvePoint::infinity();
    let approx =
        approximate_abel_jacobi_integral(&curve, &point, AbelJacobiConfig::strict()).unwrap();

    assert_eq!(approx.curve(), &curve);
    assert_eq!(approx.point(), &point);
    assert_eq!(approx.contour().radius(), 0.0);
    assert_eq!(approx.contour().tail_length(), 0.0);
    assert_eq!(approx.value(), &Complex64::new(0.0, 0.0));
    assert_eq!(
        approx.initial_branch_choice(),
        AbelJacobiInitialBranchChoice::Principal
    );
    assert_eq!(approx.segment_integral(), &Complex64::new(0.0, 0.0));
    assert_eq!(approx.ray_integral(), &Complex64::new(0.0, 0.0));
    assert_eq!(approx.tail_correction(), &Complex64::new(0.0, 0.0));
    assert_eq!(approx.integration_steps_used(), 0);
    assert_eq!(approx.branch_adjustments_used(), 0);
}

#[test]
fn approximate_abel_jacobi_integral_recovers_a_generic_square_lattice_sample() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(16).unwrap()).unwrap();
    let original_z = Complex64::new(0.2, 0.15);
    let map_result = map_torus_point_to_curve(
        &lattice,
        original_z,
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        ApproxTolerance::loose(),
    )
    .unwrap();

    let approx = approximate_abel_jacobi_integral(
        &curve,
        map_result.point(),
        AbelJacobiConfig {
            tolerance: ApproxTolerance::new(1.0e-2, 1.0e-2),
            integration_steps: 512,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 16,
            max_lattice_corrections: 4,
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            validation_policy: AbelJacobiValidationPolicy::strict(),
        },
    )
    .unwrap();

    assert!(ApproxTolerance::new(1.0e-2, 1.0e-2).real_close(approx.value().re, original_z.re));
    assert!(ApproxTolerance::new(1.0e-2, 1.0e-2).real_close(approx.value().im, original_z.im));
    assert!(matches!(
        approx.initial_branch_choice(),
        AbelJacobiInitialBranchChoice::Principal | AbelJacobiInitialBranchChoice::Alternate
    ));
    assert!((*approx.segment_integral()).is_finite());
    assert!((*approx.ray_integral()).is_finite());
    assert!((*approx.tail_correction()).is_finite());
    assert!(approx.integration_steps_used() >= 2);
}

#[test]
fn abel_jacobi_point_recovery_report_preserves_explicit_layers() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let periods = RecoveredPeriodBasis::from_lattice(ComplexLattice::from_tau(tau.clone()));
    let representative = Complex64::new(0.25, 0.5);
    let torus_point = periods
        .lattice()
        .reduce_complex_point_to_torus_point(representative)
        .unwrap();
    let contour = AbelJacobiContourReport::new(
        LegendreContourStrategy::CanonicalSegmentThenRay,
        Complex64::new(0.1, 0.2),
        Complex64::new(1.5, -0.75),
        0.7,
        2.5,
        10.0,
        0.4,
    )
    .unwrap();
    let metadata = AbelJacobiRecoveryMetadata::new(
        AbelJacobiRecoveryStatus::Succeeded,
        128,
        2,
        1,
        ApproxTolerance::strict(),
        Some(1.0e-10),
        Some(2.0e-10),
    );
    let point = AnalyticCurvePoint::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0));

    let report = AbelJacobiPointRecoveryReport::new(
        periods.clone(),
        AbelJacobiIntegralApprox::new(
            curve.clone(),
            point.clone(),
            contour.clone(),
            representative,
            AbelJacobiIntegralDecomposition::new(
                AbelJacobiInitialBranchChoice::Principal,
                Complex64::new(0.1, 0.0),
                Complex64::new(0.1, 0.0),
                Complex64::new(0.05, 0.5),
            )
            .unwrap(),
            AbelJacobiIntegralNumerics::new(128, 2, ApproxTolerance::strict()),
        )
        .unwrap(),
        torus_point,
        representative,
        AbelJacobiRoundtripValidationReport::new(
            point.clone(),
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            Some(ComplexApproxComparison::new(
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                ApproxTolerance::strict(),
            )),
            Some(ComplexApproxComparison::new(
                Complex64::new(2.0, 0.0),
                Complex64::new(2.0, 0.0),
                ApproxTolerance::strict(),
            )),
        ),
        metadata.clone(),
    )
    .unwrap();

    assert_eq!(report.curve(), &curve);
    assert_eq!(report.point(), &point);
    assert_eq!(report.periods(), &periods);
    assert_eq!(report.contour(), &contour);
    assert_eq!(report.tau(), tau);
    assert_eq!(report.raw_integral_value(), &representative);
    assert_eq!(report.reduced_representative(), &representative);
    assert_eq!(
        report.validation_report().lattice_truncation(),
        LatticeSumTruncation::new(16).unwrap()
    );
    assert_eq!(report.metadata(), &metadata);
    assert!(
        periods.lattice().torus_points_eq(
            report.torus_point(),
            &periods
                .lattice()
                .reduce_complex_point_to_torus_point(representative)
                .unwrap()
        )
    );
}

#[test]
fn point_roundtrip_validation_report_reuses_point_recovery_layers() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let periods = RecoveredPeriodBasis::from_lattice(ComplexLattice::from_tau(tau.clone()));
    let representative = Complex64::new(0.25, 0.5);
    let torus_point = periods
        .lattice()
        .reduce_complex_point_to_torus_point(representative)
        .unwrap();
    let point = AnalyticCurvePoint::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0));
    let validation_config = PointRoundTripValidationConfig::strict();
    let point_recovery_report = AbelJacobiPointRecoveryReport::new(
        periods.clone(),
        AbelJacobiIntegralApprox::new(
            curve.clone(),
            point.clone(),
            AbelJacobiContourReport::new(
                LegendreContourStrategy::CanonicalSegmentThenRay,
                Complex64::new(0.1, 0.2),
                Complex64::new(1.5, -0.75),
                0.7,
                2.5,
                10.0,
                0.4,
            )
            .unwrap(),
            representative,
            AbelJacobiIntegralDecomposition::new(
                AbelJacobiInitialBranchChoice::Principal,
                Complex64::new(0.1, 0.0),
                Complex64::new(0.1, 0.0),
                Complex64::new(0.05, 0.5),
            )
            .unwrap(),
            AbelJacobiIntegralNumerics::new(128, 2, ApproxTolerance::strict()),
        )
        .unwrap(),
        torus_point,
        representative,
        AbelJacobiRoundtripValidationReport::new(
            point.clone(),
            validation_config.lattice_truncation(),
            validation_config.function_truncation(),
            Some(ComplexApproxComparison::new(
                Complex64::new(1.0, 0.0),
                Complex64::new(1.0, 0.0),
                validation_config.tolerance(),
            )),
            Some(ComplexApproxComparison::new(
                Complex64::new(2.0, 0.0),
                Complex64::new(2.0, 0.0),
                validation_config.tolerance(),
            )),
        ),
        AbelJacobiRecoveryMetadata::new(
            AbelJacobiRecoveryStatus::Succeeded,
            128,
            2,
            0,
            ApproxTolerance::strict(),
            Some(0.0),
            Some(0.0),
        ),
    )
    .unwrap();

    let report =
        PointRoundTripValidationReport::new(point_recovery_report, validation_config).unwrap();

    assert_eq!(report.point(), &point);
    assert_eq!(report.recovered_curve_point(), &point);
    assert_eq!(
        report.lattice_truncation(),
        validation_config.lattice_truncation()
    );
    assert_eq!(
        report.function_truncation(),
        validation_config.function_truncation()
    );
    assert_eq!(report.tolerance(), validation_config.tolerance());
    assert!(report.agrees_approximately());
}

#[test]
fn inverse_uniformization_point_report_checks_curve_and_period_coherence() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(12).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let basis_report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();
    let representative = Complex64::new(0.25, 0.5);
    let torus_point = basis_report
        .periods()
        .lattice()
        .reduce_complex_point_to_torus_point(representative)
        .unwrap();
    let reduced_representative = basis_report
        .periods()
        .lattice()
        .point_from_fundamental_coordinates(torus_point.coordinate().clone());
    let metadata = AbelJacobiRecoveryMetadata::new(
        AbelJacobiRecoveryStatus::Succeeded,
        64,
        1,
        0,
        ApproxTolerance::strict(),
        None,
        None,
    );
    let point = AnalyticCurvePoint::infinity();
    let point_report = AbelJacobiPointRecoveryReport::new(
        basis_report.periods().clone(),
        AbelJacobiIntegralApprox::new(
            curve.clone(),
            point.clone(),
            AbelJacobiContourReport::new(
                LegendreContourStrategy::CanonicalSegmentThenRay,
                Complex64::new(0.25, 0.5),
                Complex64::new(0.25, 0.5),
                0.0,
                0.0,
                0.0,
                0.0,
            )
            .unwrap(),
            Complex64::new(0.0, 0.0),
            AbelJacobiIntegralDecomposition::new(
                AbelJacobiInitialBranchChoice::Principal,
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
            )
            .unwrap(),
            AbelJacobiIntegralNumerics::new(0, 0, ApproxTolerance::strict()),
        )
        .unwrap(),
        torus_point,
        reduced_representative,
        AbelJacobiRoundtripValidationReport::new(
            point.clone(),
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            None,
            None,
        ),
        metadata.clone(),
    )
    .unwrap();

    let combined =
        InverseUniformizationPointRecoveryReport::new(basis_report.clone(), point_report.clone())
            .unwrap();

    assert_eq!(combined.curve(), &curve);
    assert_eq!(combined.point(), &point);
    assert_eq!(combined.periods(), basis_report.periods());
    assert_eq!(combined.contour().radius(), 0.0);
    assert_eq!(combined.contour().tail_length(), 0.0);
    assert_eq!(combined.tau(), basis_report.tau());
    assert_eq!(combined.metadata(), &metadata);

    let wrong_periods =
        RecoveredPeriodBasis::new(Complex64::new(2.0, 0.0), Complex64::new(0.0, 2.0)).unwrap();
    let wrong_torus_point = wrong_periods
        .lattice()
        .reduce_complex_point_to_torus_point(representative)
        .unwrap();
    let wrong_point_report = AbelJacobiPointRecoveryReport::new(
        wrong_periods,
        AbelJacobiIntegralApprox::new(
            curve,
            point,
            AbelJacobiContourReport::new(
                LegendreContourStrategy::CanonicalSegmentThenRay,
                Complex64::new(0.25, 0.5),
                Complex64::new(0.25, 0.5),
                0.0,
                0.0,
                0.0,
                0.0,
            )
            .unwrap(),
            Complex64::new(0.0, 0.0),
            AbelJacobiIntegralDecomposition::new(
                AbelJacobiInitialBranchChoice::Principal,
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
            )
            .unwrap(),
            AbelJacobiIntegralNumerics::new(0, 0, ApproxTolerance::strict()),
        )
        .unwrap(),
        wrong_torus_point,
        representative,
        AbelJacobiRoundtripValidationReport::new(
            AnalyticCurvePoint::infinity(),
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            None,
            None,
        ),
        metadata,
    )
    .unwrap();

    assert!(matches!(
        InverseUniformizationPointRecoveryReport::new(basis_report, wrong_point_report),
        Err(AnalyticCurveError::InverseUniformizationFailed)
    ));
}

#[test]
fn recover_torus_point_from_curve_point_with_periods_recovers_a_generic_square_lattice_sample() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let periods = RecoveredPeriodBasis::from_lattice(lattice.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(16).unwrap()).unwrap();
    let original_z = Complex64::new(0.2, 0.15);
    let map_result = map_torus_point_to_curve(
        &lattice,
        original_z,
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        ApproxTolerance::loose(),
    )
    .unwrap();

    let recovered = recover_torus_point_from_curve_point_with_periods(
        &curve,
        map_result.point(),
        &periods,
        AbelJacobiConfig {
            tolerance: ApproxTolerance::new(1.0e-2, 1.0e-2),
            integration_steps: 512,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 16,
            max_lattice_corrections: 4,
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            validation_policy: AbelJacobiValidationPolicy::strict(),
        },
    )
    .unwrap();

    assert!(recovered.contour().radius().is_sign_positive());
    assert!(recovered.contour().tail_length().is_sign_positive());
    assert!(
        recovered
            .contour()
            .min_distance_to_branch_points()
            .is_sign_positive()
    );
    assert!(
        ApproxTolerance::new(1.0e-2, 1.0e-2)
            .real_close(recovered.reduced_representative().re, original_z.re,)
    );
    assert!(
        ApproxTolerance::new(1.0e-2, 1.0e-2)
            .real_close(recovered.reduced_representative().im, original_z.im,)
    );
    assert!(recovered.metadata().succeeded());
    assert!(recovered.validation_report().agrees_approximately());
    assert!(
        recovered
            .validation_report()
            .x_comparison()
            .unwrap()
            .agrees_approximately()
    );
    assert!(
        recovered
            .validation_report()
            .y_comparison()
            .unwrap()
            .agrees_approximately()
    );
    assert!(
        recovered
            .metadata()
            .validation_x_residual_norm()
            .unwrap()
            .is_finite()
    );
    assert!(
        recovered
            .metadata()
            .validation_y_residual_norm()
            .unwrap()
            .is_finite()
    );
}

#[test]
fn point_roundtrip_validation_with_periods_reports_successful_generic_sample() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let periods = RecoveredPeriodBasis::from_lattice(lattice.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(16).unwrap()).unwrap();
    let original_z = Complex64::new(0.2, 0.15);
    let point = map_torus_point_to_curve(
        &lattice,
        original_z,
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        ApproxTolerance::loose(),
    )
    .unwrap()
    .point()
    .clone();

    let report = validate_point_inverse_uniformization_roundtrip_with_periods(
        &curve,
        &point,
        &periods,
        AbelJacobiConfig {
            tolerance: ApproxTolerance::new(1.0e-2, 1.0e-2),
            integration_steps: 512,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 16,
            max_lattice_corrections: 4,
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            validation_policy: AbelJacobiValidationPolicy::strict(),
        },
        PointRoundTripValidationConfig::new(
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            ApproxTolerance::new(1.0e-2, 1.0e-2),
        ),
    )
    .unwrap();

    assert!(report.agrees_approximately());
    assert!(
        ApproxTolerance::new(1.0e-2, 1.0e-2)
            .real_close(report.reduced_representative().re, original_z.re)
    );
    assert!(
        ApproxTolerance::new(1.0e-2, 1.0e-2)
            .real_close(report.reduced_representative().im, original_z.im)
    );
    assert!(report.x_comparison().unwrap().agrees_approximately());
    assert!(report.y_comparison().unwrap().agrees_approximately());
}

#[test]
fn point_roundtrip_validation_can_return_a_failed_report_instead_of_erroring() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let periods = RecoveredPeriodBasis::from_lattice(lattice.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(18).unwrap()).unwrap();
    let point = map_torus_point_to_curve(
        &lattice,
        Complex64::new(0.23, 0.17),
        LatticeSumTruncation::new(18).unwrap(),
        EllipticFunctionTruncation::new(16).unwrap(),
        ApproxTolerance::strict(),
    )
    .unwrap()
    .point()
    .clone();

    let validation_config = PointRoundTripValidationConfig::new(
        LatticeSumTruncation::new(1).unwrap(),
        EllipticFunctionTruncation::new(1).unwrap(),
        ApproxTolerance::strict(),
    );

    let report = validate_point_inverse_uniformization_roundtrip_with_periods(
        &curve,
        &point,
        &periods,
        AbelJacobiConfig {
            tolerance: ApproxTolerance::strict(),
            validation_policy: AbelJacobiValidationPolicy::strict(),
            ..AbelJacobiConfig::strict()
        },
        validation_config,
    )
    .unwrap();

    assert!(!report.agrees_approximately());
    assert_eq!(
        report.point_recovery_report().metadata().status(),
        AbelJacobiRecoveryStatus::ValidationFailed
    );
}

#[test]
fn approximate_abel_jacobi_integral_rejects_branch_points_with_y_near_zero() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let roots = recover_weierstrass_cubic_roots(&curve, PeriodRecoveryConfig::strict()).unwrap();
    let point = AnalyticCurvePoint::new(*roots.roots()[0], Complex64::new(0.0, 0.0));

    assert!(matches!(
        approximate_abel_jacobi_integral(&curve, &point, AbelJacobiConfig::strict()),
        Err(AnalyticCurveError::AbelJacobiIntegrationFailed)
    ));
}

#[test]
fn roundtrip_validation_can_fail_when_validation_budget_is_too_small() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let periods = RecoveredPeriodBasis::from_lattice(lattice.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(18).unwrap()).unwrap();
    let point = map_torus_point_to_curve(
        &lattice,
        Complex64::new(0.23, 0.17),
        LatticeSumTruncation::new(18).unwrap(),
        EllipticFunctionTruncation::new(16).unwrap(),
        ApproxTolerance::strict(),
    )
    .unwrap()
    .point()
    .clone();

    let config = AbelJacobiConfig {
        tolerance: ApproxTolerance::strict(),
        validation_policy: AbelJacobiValidationPolicy {
            lattice_truncation_radius: 1,
            function_truncation_radius: 1,
        },
        ..AbelJacobiConfig::strict()
    };

    assert!(matches!(
        recover_torus_point_from_curve_point_with_periods(&curve, &point, &periods, config),
        Err(AnalyticCurveError::PeriodValidationFailed)
    ));
}

#[test]
fn end_to_end_inverse_uniformization_entrypoint_recovers_the_identity_class() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let point = AnalyticCurvePoint::infinity();
    let report = recover_torus_point_from_curve_point(
        &curve,
        &point,
        PeriodRecoveryConfig::loose(),
        AbelJacobiConfig::loose(),
    )
    .unwrap();

    assert_eq!(report.point(), &point);
    assert_eq!(report.reduced_representative(), &Complex64::new(0.0, 0.0));
    assert_eq!(report.metadata().validation_x_residual_norm(), Some(0.0));
    assert_eq!(report.metadata().validation_y_residual_norm(), Some(0.0));
}

#[test]
fn end_to_end_point_roundtrip_validation_handles_the_identity_class() {
    let tau = UpperHalfPlanePoint::tau_i();
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap()).unwrap();
    let point = AnalyticCurvePoint::infinity();
    let report = validate_point_inverse_uniformization_roundtrip(
        &curve,
        &point,
        PeriodRecoveryConfig::loose(),
        AbelJacobiConfig::loose(),
        PointRoundTripValidationConfig::loose(),
    )
    .unwrap();

    assert_eq!(report.point(), &point);
    assert_eq!(report.recovered_curve_point(), &point);
    assert_eq!(report.reduced_representative(), &Complex64::new(0.0, 0.0));
    assert!(report.agrees_approximately());
}

#[test]
fn end_to_end_point_roundtrip_with_recovered_periods_succeeds_for_a_square_lattice_sample() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau.clone());
    let curve =
        AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(18).unwrap()).unwrap();
    let point = map_torus_point_to_curve(
        &lattice,
        Complex64::new(0.2, 0.15),
        LatticeSumTruncation::new(18).unwrap(),
        EllipticFunctionTruncation::new(16).unwrap(),
        ApproxTolerance::strict(),
    )
    .unwrap()
    .point()
    .clone();

    let report = validate_point_inverse_uniformization_roundtrip(
        &curve,
        &point,
        PeriodRecoveryConfig::strict(),
        AbelJacobiConfig {
            tolerance: ApproxTolerance::new(1.0e-2, 1.0e-2),
            integration_steps: 512,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 16,
            max_lattice_corrections: 4,
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            validation_policy: AbelJacobiValidationPolicy::strict(),
        },
        PointRoundTripValidationConfig::new(
            LatticeSumTruncation::new(16).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            ApproxTolerance::new(1.0e-2, 1.0e-2),
        ),
    )
    .unwrap();

    assert!(report.agrees_approximately());
    assert_eq!(
        report.point_recovery_report().metadata().status(),
        AbelJacobiRecoveryStatus::Succeeded
    );
}

#[test]
fn inverse_uniformization_rejects_point_not_on_curve() {
    let curve =
        AnalyticWeierstrassCurve::new(Complex64::new(0.0, 0.0), Complex64::new(4.0, 0.0)).unwrap();
    let off_curve_point =
        AnalyticCurvePoint::new(Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0));

    assert!(!curve.contains_point_approx(&off_curve_point, ApproxTolerance::strict()));
    assert!(matches!(
        recover_torus_point_from_curve_point(
            &curve,
            &off_curve_point,
            PeriodRecoveryConfig::strict(),
            AbelJacobiConfig::strict(),
        ),
        Err(AnalyticCurveError::InverseUniformizationFailed)
            | Err(AnalyticCurveError::PeriodValidationFailed)
            | Err(AnalyticCurveError::AbelJacobiIntegrationFailed)
            | Err(AnalyticCurveError::BranchChoiceAmbiguous)
    ));
}
