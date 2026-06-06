    use num_complex::Complex64;

    use crate::visualization::elliptic_curves::analytic::{
        describe_analytic_division_polynomial_comparison,
        describe_analytic_even_division_polynomial_report, describe_analytic_invariants,
        describe_analytic_odd_division_polynomial_report, describe_analytic_torsion_point_approx,
        describe_canonical_tau_recovery_report, describe_complex_lattice,
        describe_cubic_root_configuration_report, describe_cubic_root_recovery_report,
        describe_eisenstein_sum, describe_fundamental_domain_reduction_report,
        describe_fundamental_domain_reduction_step, describe_invariant_recovery_validation_report,
        describe_inverse_uniformization_j_validation_report, describe_j_invariant_comparison,
        describe_legendre_parameter, describe_legendre_parameter_conditioning,
        describe_legendre_parameter_orbit, describe_legendre_reduction,
        describe_legendre_reduction_report, describe_modular_invariance_report,
        describe_modular_matrix, describe_numerical_recovery_metadata,
        describe_period_basis_recovery_report, describe_period_lattice,
        describe_period_recovery_config, describe_period_recovery_report,
        describe_point_roundtrip_validation_config, describe_point_roundtrip_validation_report,
        describe_q_parameter, describe_recovered_period_basis,
        describe_recovered_period_basis_report, describe_tau_recovery_report,
        describe_torus_to_curve_map, describe_weierstrass_cubic_roots,
        describe_weierstrass_differential_equation, describe_weierstrass_p_approx,
        format_analytic_cubic_model, format_complex_scalar_compact,
        format_short_weierstrass_over_complex,
    };
    use crate::elliptic_curves::analytic::{
        AbelJacobiConfig, AbelJacobiValidationPolicy, LegendreContourStrategy,
        PointRoundTripValidationConfig,
        validate_point_inverse_uniformization_roundtrip_with_periods,
    };
    use crate::elliptic_curves::{
        AnalyticCurvePoint, AnalyticDivisionPolynomialComparisonCase, AnalyticWeierstrassCurve,
        ApproxTolerance, ComplexLattice, EllipticFunctionTruncation, LatticeSumTruncation,
        LegendreParameter, LegendreParameterConditioning, LegendreReduction,
        LegendreReductionReport, ModularMatrix, ModularQParameter, NumericalRecoveryMetadata,
        PeriodLatticeApprox, PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryStatus,
        QExpansionTruncation, RecoveredPeriodBasis, UpperHalfPlanePoint, WeierstrassCubicRoots,
        analytic_invariants, compare_analytic_torsion_with_division_polynomial,
        compare_j_from_eisenstein_and_q_expansion,
        compare_primitive_analytic_torsion_with_division_polynomial,
        cubic_root_configuration_report, g4_sum, map_torus_point_to_curve,
        recover_canonical_tau_from_curve, recover_period_basis, recover_tau_from_curve,
        recover_weierstrass_cubic_roots_with_report, reduce_tau_to_standard_fundamental_domain,
        validate_recovered_lattice_invariants, validate_recovered_tau_by_j_invariant,
        verify_j_modular_invariance, verify_weierstrass_differential_equation, weierstrass_p,
    };
    use crate::visualization::Visualizable;
    use crate::visualization::elliptic_curves::format_point_compact;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn lattice_description_mentions_basis_and_tau() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_complex_lattice(&lattice);

        assert!(text.contains("Complex lattice"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
    }

    #[test]
    fn eisenstein_description_mentions_weight_and_truncation() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let sum = g4_sum(&lattice, LatticeSumTruncation::default_educational()).unwrap();
        let text = describe_eisenstein_sum(&sum);

        assert!(text.contains("Eisenstein sum"));
        assert!(text.contains("weight k = 4"));
        assert!(text.contains("truncation radius"));
        assert!(text.contains("value"));
    }

    #[test]
    fn analytic_invariant_description_mentions_g2_g3_delta_and_j() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let invariants =
            analytic_invariants(&lattice, LatticeSumTruncation::default_educational()).unwrap();
        let text = describe_analytic_invariants(&invariants);

        assert!(text.contains("Analytic invariants"));
        assert!(text.contains("g₂"));
        assert!(text.contains("g₃"));
        assert!(text.contains("Δ"));
        assert!(text.contains("j"));
    }

    #[test]
    fn q_parameter_description_mentions_tau_q_and_open_unit_disc() {
        let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_q_parameter(&q);

        assert!(text.contains("Modular q-parameter"));
        assert!(text.contains("τ ="));
        assert!(text.contains("q = e^(2πiτ)"));
        assert!(text.contains("|q|"));
        assert!(text.contains("open unit disc"));
    }

    #[test]
    fn period_recovery_config_description_mentions_all_budgets() {
        let config = PeriodRecoveryConfig::educational_default();
        let text = describe_period_recovery_config(&config);

        assert!(text.contains("Period recovery config"));
        assert!(text.contains("Newton iteration budget"));
        assert!(text.contains("AGM iteration budget"));
        assert!(text.contains("Abel-Jacobi integration steps"));
        assert!(text.contains("branch lattice search radius"));
        assert!(text.contains("fundamental-domain reduction steps"));
    }

    #[test]
    fn period_lattice_description_mentions_basis_and_tau() {
        let periods = PeriodLatticeApprox::standard_from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_period_lattice(&periods);

        assert!(text.contains("Approximate period lattice"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
        assert!(text.contains("not a canonical lattice representative"));
    }

    #[test]
    fn numerical_recovery_metadata_description_mentions_method_status_and_counters() {
        let metadata = NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::ValidationFailed,
            7,
            0,
            0,
            2,
            ApproxTolerance::strict(),
            Some(1.0e-9),
        )
        .with_cardano_diagnostics(Complex64::new(3.0, -4.0), 2.5e-14, 0, 2);
        let text = describe_numerical_recovery_metadata(&metadata);

        assert!(text.contains("Numerical recovery metadata"));
        assert!(text.contains("resolved method = hybrid"));
        assert!(text.contains("status = validation failed"));
        assert!(text.contains("newton iterations used = 7"));
        assert!(text.contains("validation residual norm"));
        assert!(text.contains("Cardano discriminant"));
        assert!(text.contains("Cardano branch product residual norm"));
        assert!(text.contains("selected Cardano branch indices"));
        assert!(text.contains("used principal Cardano branches = no"));
    }

    #[test]
    fn legendre_parameter_description_mentions_lambda_and_orbit_caveat() {
        let parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let text = describe_legendre_parameter(&parameter);

        assert!(text.contains("Legendre parameter"));
        assert!(text.contains("lambda"));
        assert!(text.contains("1 - lambda"));
        assert!(text.contains("six-element S3 orbit"));
    }

    #[test]
    fn legendre_parameter_description_avoids_compact_rounding_near_zero() {
        let parameter = LegendreParameter::new(c(-3.333333333e-8, 0.0)).unwrap();
        let text = describe_legendre_parameter(&parameter);

        assert!(text.contains("lambda ≈ -3.333333333"));
        assert!(text.contains("e-8"));
        assert!(!text.contains("lambda ≈ 0\n"));
    }

    #[test]
    fn legendre_parameter_orbit_description_lists_all_classical_transforms() {
        let parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let text = describe_legendre_parameter_orbit(&parameter.orbit());

        assert!(text.contains("Legendre parameter orbit"));
        assert!(text.contains("lambda ≈"));
        assert!(text.contains("1 - lambda"));
        assert!(text.contains("1 / lambda"));
        assert!(text.contains("same Legendre class"));
    }

    #[test]
    fn legendre_conditioning_description_mentions_singularity_verdict() {
        let text =
            describe_legendre_parameter_conditioning(LegendreParameterConditioning::NearInfinity);

        assert!(text.contains("Legendre parameter conditioning"));
        assert!(text.contains("near infinity"));
        assert!(text.contains("near singular locus = yes"));
    }

    #[test]
    fn legendre_reduction_description_mentions_permutation_and_principal_branch_scales() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let text = describe_legendre_reduction(&reduction);

        assert!(text.contains("Legendre reduction"));
        assert!(text.contains("selected permutation"));
        assert!(text.contains("not a canonical root ordering"));
        assert!(text.contains("rhs scale factor"));
        assert!(text.contains("principal sqrt(x scale)"));
        assert!(text.contains("principal y scale"));
        assert!(text.contains("invariant differential scale"));
    }

    #[test]
    fn legendre_reduction_report_description_mentions_input_order_caveat() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let report =
            LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let text = describe_legendre_reduction_report(&report);

        assert!(text.contains("Legendre reduction report"));
        assert!(text.contains("selected orbit element relative to input order"));
        assert!(text.contains("conditioning = generic"));
        assert!(text.contains("singularity distance"));
        assert!(text.contains("maximizing min(|lambda|, |1 - lambda|, 1 / |lambda|)"));
        assert!(text.contains("not canonical by itself"));
        assert!(text.contains("reduction summary"));
    }

    #[test]
    fn weierstrass_cubic_roots_description_mentions_roots_and_invariants() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_weierstrass_cubic_roots(&roots);

        assert!(text.contains("Weierstrass cubic roots"));
        assert!(text.contains("root[0]"));
        assert!(text.contains("not canonical"));
        assert!(text.contains("g₂"));
        assert!(text.contains("g₃"));
        assert!(text.contains("minimum pairwise distance"));
    }

    #[test]
    fn weierstrass_cubic_roots_description_uses_diagnostic_precision_for_nearly_colliding_roots() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(1.0 + 1.0e-7, 0.0),
            c(-2.0 - 1.0e-7, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_weierstrass_cubic_roots(&roots);

        assert!(text.contains("root[0] ≈ 1"));
        assert!(text.contains("root[1] ≈ 1.0000001"));
        assert!(text.contains("root[2] ≈ -2.0000001"));
    }

    #[test]
    fn cubic_root_configuration_description_mentions_shape_and_separation() {
        let roots = WeierstrassCubicRoots::new(
            c(2.0, 1.0),
            c(-3.0, 0.0),
            c(2.0, -1.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let report = cubic_root_configuration_report(&roots, ApproxTolerance::strict());
        let text = describe_cubic_root_configuration_report(&report);

        assert!(text.contains("Cubic-root configuration"));
        assert!(text.contains("approximately conjugate pair"));
        assert!(text.contains("separation = well separated"));
        assert!(text.contains("best conjugate-pair residual"));
        assert!(text.contains("roots summary"));
    }

    #[test]
    fn j_invariant_comparison_description_mentions_both_routes_and_difference() {
        let report = compare_j_from_eisenstein_and_q_expansion(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::default_educational(),
            QExpansionTruncation::new(3).unwrap(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_j_invariant_comparison(&report);

        assert!(text.contains("j-invariant comparison"));
        assert!(text.contains("j from Eisenstein sums"));
        assert!(text.contains("j from q-expansion"));
        assert!(text.contains("|difference|"));
        assert!(text.contains("agrees under tolerance"));
    }

    #[test]
    fn period_recovery_report_description_mentions_periods_and_j_residual() {
        let tau = UpperHalfPlanePoint::tau_i();
        let curve =
            AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap())
                .unwrap();
        let periods = PeriodLatticeApprox::standard_from_tau(tau);
        let recovered_j = curve.j_invariant().unwrap();
        let report = crate::elliptic_curves::PeriodRecoveryReport::new(
            curve,
            periods,
            recovered_j,
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_period_recovery_report(&report);

        assert!(text.contains("Period recovery report"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("recovered j"));
        assert!(text.contains("curve-side j"));
        assert!(text.contains("|difference|"));
    }

    #[test]
    fn inverse_uniformization_j_validation_description_mentions_tau_invariants_and_residual() {
        let tau = UpperHalfPlanePoint::tau_i();
        let truncation = LatticeSumTruncation::new(12).unwrap();
        let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
        let report = validate_recovered_tau_by_j_invariant(
            &curve,
            &tau,
            truncation,
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_inverse_uniformization_j_validation_report(&report);

        assert!(text.contains("Inverse-uniformization j-validation report"));
        assert!(text.contains("τ ≈"));
        assert!(text.contains("recovered g₂"));
        assert!(text.contains("recovered g₃"));
        assert!(text.contains("recovered j"));
        assert!(text.contains("curve-side j"));
        assert!(text.contains("lattice truncation radius"));
        assert!(text.contains("modular j-class"));
    }

    #[test]
    fn invariant_recovery_validation_description_mentions_interpretation_and_all_invariants() {
        let tau = UpperHalfPlanePoint::tau_i();
        let truncation = LatticeSumTruncation::new(12).unwrap();
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
        let text = describe_invariant_recovery_validation_report(&report);

        assert!(text.contains("Invariant recovery validation report"));
        assert!(text.contains("recovered g₂"));
        assert!(text.contains("curve-side g₂"));
        assert!(text.contains("recovered g₃"));
        assert!(text.contains("curve-side g₃"));
        assert!(text.contains("recovered Δ"));
        assert!(text.contains("curve-side Δ"));
        assert!(text.contains("interpretation ="));
        assert!(text.contains("same modular class via j"));
        assert!(text.contains("homothety-invariant"));
    }

    #[test]
    fn point_roundtrip_validation_config_description_mentions_both_truncations() {
        let config = PointRoundTripValidationConfig::strict();
        let text = describe_point_roundtrip_validation_config(&config);

        assert!(text.contains("Point roundtrip validation config"));
        assert!(text.contains("lattice truncation radius"));
        assert!(text.contains("elliptic-function truncation radius"));
        assert!(text.contains("tolerance ="));
        assert!(text.contains("forward check z -> (wp(z), wp'(z))"));
    }

    #[test]
    fn point_roundtrip_validation_report_description_mentions_torus_and_curve_sides() {
        let tau = UpperHalfPlanePoint::tau_i();
        let lattice = ComplexLattice::from_tau(tau.clone());
        let periods = RecoveredPeriodBasis::from_lattice(lattice.clone());
        let curve =
            AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(16).unwrap())
                .unwrap();
        let point = map_torus_point_to_curve(
            &lattice,
            Complex64::new(0.2, 0.15),
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
        let text = describe_point_roundtrip_validation_report(&report);

        assert!(text.contains("Point inverse-uniformization roundtrip report"));
        assert!(text.contains("source point P ="));
        assert!(text.contains("recovered torus representative"));
        assert!(text.contains("torus coordinates in the recovered basis"));
        assert!(text.contains("recovered curve point"));
        assert!(text.contains("x residual norm"));
        assert!(text.contains("y residual norm"));
        assert!(text.contains("lattice truncation radius"));
        assert!(text.contains("elliptic-function truncation radius"));
        assert!(text.contains("reuses the recovered torus class"));
    }

    #[test]
    fn recovered_period_basis_description_mentions_basis_tau_and_covolume() {
        let basis =
            RecoveredPeriodBasis::new(c(2.0, 0.0), c(1.0, 3.0)).expect("valid positive basis");
        let text = describe_recovered_period_basis(&basis);

        assert!(text.contains("Recovered period basis"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
        assert!(text.contains("covolume"));
        assert!(text.contains("not a canonical"));
    }

    #[test]
    fn recovered_period_basis_report_description_mentions_transport_story() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let period_basis_report = recover_period_basis(
            &AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap(),
            PeriodRecoveryConfig::strict(),
        )
        .unwrap();
        let text = describe_recovered_period_basis_report(period_basis_report.basis_report());

        assert!(text.contains("Recovered period basis report"));
        assert!(text.contains("basis summary"));
        assert!(text.contains("invariant differential scale"));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("complete elliptic integral summary"));
        assert!(text.contains("transporting the normalized Legendre periods"));
        assert!(text.contains("τ = ω₂ / ω₁"));
    }

    #[test]
    fn period_basis_recovery_report_description_mentions_curve_roots_tau_and_metadata() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_period_basis_recovery_report(&report);

        assert!(text.contains("Period-basis recovery report"));
        assert!(text.contains("curve ="));
        assert!(text.contains("roots ="));
        assert!(text.contains("configuration ="));
        assert!(text.contains("separation ="));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("period basis summary"));
        assert!(text.contains("τ summary"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("not canonical"));
    }

    #[test]
    fn tau_recovery_report_description_mentions_tau_and_no_parallel_pipeline_caveat() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report = recover_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_tau_recovery_report(&report);

        assert!(text.contains("Tau recovery report"));
        assert!(text.contains("τ ≈"));
        assert!(text.contains("period basis summary"));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("not a separate tau-only algorithm"));
    }

    #[test]
    fn canonical_tau_recovery_report_description_mentions_original_tau_canonical_tau_and_matrix() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report =
            recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_canonical_tau_recovery_report(&report);

        assert!(text.contains("Canonical tau recovery report"));
        assert!(text.contains("original τ"));
        assert!(text.contains("canonical τ"));
        assert!(text.contains("accumulated modular matrix"));
        assert!(text.contains("fundamental-domain status"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("applying the accumulated modular matrix"));
    }

    #[test]
    fn cubic_root_recovery_report_description_mentions_reconstruction_and_metadata() {
        let curve = AnalyticWeierstrassCurve::new(c(28.0, 0.0), c(-24.0, 0.0)).unwrap();
        let report =
            recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
                .unwrap();
        let text = describe_cubic_root_recovery_report(&report);

        assert!(text.contains("Cubic-root recovery report"));
        assert!(text.contains("configuration ="));
        assert!(text.contains("separation ="));
        assert!(text.contains("reconstructed g₂"));
        assert!(text.contains("curve-side g₃"));
        assert!(text.contains("metadata summary"));
    }

    #[test]
    fn modular_matrix_description_mentions_entries_and_action() {
        let text = describe_modular_matrix(&ModularMatrix::s());

        assert!(text.contains("Modular matrix"));
        assert!(text.contains("γ = [[0, -1], [1, 0]]"));
        assert!(text.contains("determinant = 1"));
        assert!(text.contains("action on τ"));
    }

    #[test]
    fn modular_invariance_description_mentions_both_taus_and_difference() {
        let report = verify_j_modular_invariance(
            UpperHalfPlanePoint::tau_i(),
            ModularMatrix::s(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_modular_invariance_report(&report);

        assert!(text.contains("Modular invariance check"));
        assert!(text.contains("original τ ="));
        assert!(text.contains("transformed τ ="));
        assert!(text.contains("j(τ)"));
        assert!(text.contains("j(γτ)"));
        assert!(text.contains("|difference|"));
    }

    #[test]
    fn fundamental_domain_descriptions_mention_status_and_reason() {
        let report = reduce_tau_to_standard_fundamental_domain(
            UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap(),
            8,
        )
        .unwrap();
        let report_text = describe_fundamental_domain_reduction_report(&report);
        let step_text = describe_fundamental_domain_reduction_step(&report.steps()[0]);

        assert!(report_text.contains("Fundamental-domain reduction"));
        assert!(report_text.contains("status = reduced"));
        assert!(report_text.contains("steps used ="));
        assert!(step_text.contains("Fundamental-domain reduction step"));
        assert!(step_text.contains("reason = real part lay outside the centered strip"));
    }

    #[test]
    fn weierstrass_p_description_mentions_pole_distance() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let approximation = weierstrass_p(
            &lattice,
            c(0.2, 0.15),
            EllipticFunctionTruncation::default_educational(),
        )
        .unwrap();
        let text = describe_weierstrass_p_approx(&approximation);

        assert!(text.contains("Weierstrass"));
        assert!(text.contains("truncation radius"));
        assert!(text.contains("nearest inspected pole distance"));
    }

    #[test]
    fn torus_to_curve_description_distinguishes_the_pole_case() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let map = map_torus_point_to_curve(
            &lattice,
            c(0.0, 0.0),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_torus_to_curve_map(&map);

        assert!(text.contains("Torus to curve map"));
        assert!(text.contains("values = Pole"));
        assert!(text.contains("infinity"));
    }

    #[test]
    fn analytic_torsion_point_description_mentions_index_z_and_curve_point() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let mapped = crate::elliptic_curves::map_torus_torsion_to_curve(
            &lattice,
            3,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_analytic_torsion_point_approx(&mapped[1]);

        assert!(text.contains("Analytic torsion point"));
        assert!(text.contains("torus index ="));
        assert!(text.contains("z ="));
        assert!(text.contains("curve point ="));
    }

    #[test]
    fn analytic_division_polynomial_description_distinguishes_the_pole_case() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_analytic_division_polynomial_comparison(&reports[0]);

        assert!(text.contains("Analytic torsion vs division polynomial"));
        assert!(text.contains("case = pole at identity"));
        assert!(text.contains("no finite x = ℘(z) value is available"));
    }

    #[test]
    fn odd_division_polynomial_description_mentions_psi_n_and_status() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_primitive_analytic_torsion_with_division_polynomial(
            &lattice,
            3,
            LatticeSumTruncation::larger_for_comparison(),
            EllipticFunctionTruncation::new(6).unwrap(),
            ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        let odd_report = match &reports[0] {
            AnalyticDivisionPolynomialComparisonCase::Odd(odd_report) => odd_report,
            other => panic!("expected odd report, got {other:?}"),
        };
        let text = describe_analytic_odd_division_polynomial_report(odd_report);

        assert!(text.contains("odd n"));
        assert!(text.contains("ψ_n(x)"));
        assert!(text.contains("status ="));
    }

    #[test]
    fn even_division_polynomial_description_mentions_branch_and_epsilon_n() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_primitive_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();

        let even_report = match &reports[0] {
            AnalyticDivisionPolynomialComparisonCase::Even(even_report) => even_report,
            other => panic!("expected even report, got {other:?}"),
        };
        let text = describe_analytic_even_division_polynomial_report(even_report);

        assert!(text.contains("even n"));
        assert!(text.contains("ε_n(x)"));
        assert!(text.contains("branch ="));
        assert!(text.contains("neither y(P) nor ε_n(x(P)) is approximately zero"));
    }

    #[test]
    fn differential_equation_description_mentions_lhs_rhs_and_status() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_weierstrass_differential_equation(&report);

        assert!(text.contains("Weierstrass differential equation"));
        assert!(text.contains("lhs"));
        assert!(text.contains("rhs"));
        assert!(text.contains("status"));
    }

    #[test]
    fn visualizable_trait_is_hooked_up_for_analytic_reports() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_i());
        let period_config = PeriodRecoveryConfig::educational_default();
        let period_lattice = PeriodLatticeApprox::standard_from_tau(UpperHalfPlanePoint::tau_i());
        let metadata = NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::Succeeded,
            3,
            0,
            0,
            0,
            ApproxTolerance::strict(),
            Some(1.0e-12),
        );
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let legendre_parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let legendre_orbit = legendre_parameter.orbit();
        let legendre_reduction =
            LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let legendre_report =
            LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let root_configuration = cubic_root_configuration_report(&roots, ApproxTolerance::strict());
        let map = map_torus_point_to_curve(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let point_roundtrip_config = PointRoundTripValidationConfig::loose();

        assert!(lattice.format_compact().contains("Λ = ℤ"));
        assert!(q.format_compact().contains("q(τ)"));
        assert!(period_config.describe().contains("Period recovery config"));
        assert!(
            period_lattice
                .describe()
                .contains("Approximate period lattice")
        );
        assert!(metadata.describe().contains("Numerical recovery metadata"));
        assert!(legendre_parameter.describe().contains("Legendre parameter"));
        assert!(
            legendre_orbit
                .describe()
                .contains("Legendre parameter orbit")
        );
        assert!(
            LegendreParameterConditioning::NearInfinity
                .describe()
                .contains("Legendre parameter conditioning")
        );
        assert!(legendre_reduction.describe().contains("Legendre reduction"));
        assert!(
            legendre_report
                .describe()
                .contains("Legendre reduction report")
        );
        assert!(roots.describe().contains("Weierstrass cubic roots"));
        assert!(
            root_configuration
                .describe()
                .contains("Cubic-root configuration")
        );
        assert!(q.describe().contains("Modular q-parameter"));
        assert!(ModularMatrix::s().describe().contains("Modular matrix"));
        assert!(map.describe().contains("Torus to curve map"));
        assert!(
            report
                .describe()
                .contains("Weierstrass differential equation")
        );
        assert!(
            point_roundtrip_config
                .describe()
                .contains("Point roundtrip validation config")
        );
        let modular_report = verify_j_modular_invariance(
            UpperHalfPlanePoint::tau_i(),
            ModularMatrix::s(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        assert!(
            modular_report
                .describe()
                .contains("Modular invariance check")
        );
        let reduction = reduce_tau_to_standard_fundamental_domain(
            UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap(),
            8,
        )
        .unwrap();
        assert!(
            reduction
                .describe()
                .contains("Fundamental-domain reduction")
        );
        let torsion_comparison = compare_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        assert!(
            torsion_comparison[0]
                .describe()
                .contains("Analytic torsion vs division polynomial")
        );
        let curve = AnalyticWeierstrassCurve::new(c(28.0, 0.0), c(-24.0, 0.0)).unwrap();
        let cubic_root_report =
            recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
                .unwrap();
        assert!(
            cubic_root_report
                .describe()
                .contains("Cubic-root recovery report")
        );
        let point = map_torus_point_to_curve(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap()
        .point()
        .clone();
        let point_roundtrip_report = validate_point_inverse_uniformization_roundtrip_with_periods(
            &AnalyticWeierstrassCurve::from_tau(
                &UpperHalfPlanePoint::tau_i(),
                LatticeSumTruncation::default_educational(),
            )
            .unwrap(),
            &point,
            &RecoveredPeriodBasis::from_lattice(lattice.clone()),
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
                LatticeSumTruncation::default_educational(),
                EllipticFunctionTruncation::default_educational(),
                ApproxTolerance::new(1.0e-2, 1.0e-2),
            ),
        )
        .unwrap();
        assert!(
            point_roundtrip_report
                .describe()
                .contains("Point inverse-uniformization roundtrip report")
        );
        let infinity = AnalyticCurvePoint::infinity();
        assert_eq!(format_point_compact(&infinity), "O");
    }

    #[test]
    fn specialized_curve_formatters_drop_near_zero_terms_and_imaginary_noise() {
        let analytic =
            AnalyticWeierstrassCurve::new(c(188.94472, -1.0e-15), c(1.0e-15, 2.0e-16)).unwrap();
        let short = analytic.as_short_weierstrass();

        assert_eq!(
            format_analytic_cubic_model(&analytic),
            "y^2 = 4x^3 - 188.944720x"
        );
        assert_eq!(
            format_short_weierstrass_over_complex(&short),
            "y^2 = x^3 - 47.236180x"
        );
    }

    #[test]
    fn compact_complex_formatter_drops_tiny_real_noise_next_to_large_imaginary_part() {
        let value = c(5.0e-7, 60690.762066);

        assert_eq!(format_complex_scalar_compact(&value), "60690.762066i");
        assert_eq!(format_complex_scalar_compact(&c(0.0, 0.0)), "0");
    }
