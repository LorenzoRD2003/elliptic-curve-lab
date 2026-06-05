use num_complex::Complex64;

use super::{
    CubicRootConfiguration, CubicRootConfigurationReport, CubicRootRecoveryReport,
    CubicRootSeparation, LegendreOrbitElementKind, LegendreParameter,
    LegendreParameterConditioning, LegendreReduction, LegendreReductionReport,
    NumericalRecoveryMetadata, PeriodLatticeApprox, PeriodRecoveryConfig, PeriodRecoveryMethod,
    PeriodRecoveryReport, PeriodRecoveryStatus, WeierstrassCubicRoots,
    classify_cubic_root_configuration, classify_legendre_parameter_conditioning,
    cubic_root_configuration_report, legendre_reduction_report, recover_weierstrass_cubic_roots,
    recover_weierstrass_cubic_roots_from_invariants, recover_weierstrass_cubic_roots_with_report,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
    HasAnalyticLatticeContext, HasComplexApproxComparison, LatticeSumTruncation,
    UpperHalfPlanePoint,
};
use crate::fields::ComplexApprox;

#[test]
fn config_constructor_preserves_caller_supplied_values() {
    let tolerance = ApproxTolerance::new(1.0e-8, 2.0e-8);
    let config = PeriodRecoveryConfig::new(tolerance, 9, 7, 192, 3).unwrap();

    assert_eq!(config.tolerance(), tolerance);
    assert_eq!(config.newton_max_iterations(), 9);
    assert_eq!(config.agm_max_iterations(), 7);
    assert_eq!(config.abel_jacobi_integration_steps(), 192);
    assert_eq!(config.branch_lattice_search_radius(), 3);
}

#[test]
fn config_rejects_zero_budgets() {
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 0, 1, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 0, 1, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 0, 1),
        Err(AnalyticCurveError::InvalidPeriodRecoveryConfig)
    );
    assert_eq!(
        PeriodRecoveryConfig::new(ApproxTolerance::strict(), 1, 1, 1, 0),
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

    assert_eq!(strict.tolerance(), ApproxTolerance::strict());
    assert!(strict.newton_max_iterations() > educational.newton_max_iterations());
    assert!(strict.agm_max_iterations() > educational.agm_max_iterations());
    assert!(strict.abel_jacobi_integration_steps() > educational.abel_jacobi_integration_steps());
    assert!(strict.branch_lattice_search_radius() > educational.branch_lattice_search_radius());

    assert_eq!(loose.tolerance(), ApproxTolerance::loose());
    assert!(loose.newton_max_iterations() < educational.newton_max_iterations());
    assert!(loose.agm_max_iterations() < educational.agm_max_iterations());
    assert!(loose.abel_jacobi_integration_steps() < educational.abel_jacobi_integration_steps());
    assert!(loose.branch_lattice_search_radius() < educational.branch_lattice_search_radius());
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
