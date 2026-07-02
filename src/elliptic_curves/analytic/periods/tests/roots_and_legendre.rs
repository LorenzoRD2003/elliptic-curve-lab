use super::*;
use num_complex::Complex64;

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
        LegendreParameter::new(Complex64::new(-0.25, 0.0))
            .unwrap()
            .conditioning(tolerance),
        LegendreParameterConditioning::Generic
    );
    assert_eq!(
        LegendreParameter::new(Complex64::new(1.0e-13, 0.0))
            .unwrap()
            .conditioning(tolerance),
        LegendreParameterConditioning::NearZero
    );
    assert_eq!(
        LegendreParameter::new(Complex64::new(1.0 + 1.0e-13, 0.0))
            .unwrap()
            .conditioning(tolerance),
        LegendreParameterConditioning::NearOne
    );
    assert_eq!(
        LegendreParameter::new(Complex64::new(1.0e13, 0.0))
            .unwrap()
            .conditioning(tolerance),
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
        roots
            .legendre_reduction_report(ApproxTolerance::strict())
            .unwrap(),
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
    let report = roots.configuration_report(ApproxTolerance::strict());

    assert_eq!(
        roots.classify_configuration(ApproxTolerance::strict()),
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
    let report = roots.configuration_report(ApproxTolerance::strict());

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
    let report = roots.configuration_report(ApproxTolerance::loose());

    assert_eq!(
        report.configuration(),
        CubicRootConfiguration::ThreeApproximatelyReal
    );
    assert_eq!(report.separation(), CubicRootSeparation::NearlyRepeated);
    assert!(ApproxTolerance::strict().real_close(report.min_pairwise_distance(), 5.0e-8));
}
