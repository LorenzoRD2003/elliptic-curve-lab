use super::*;

#[test]
fn recovering_periods_from_a_curve_built_from_tau_recovers_the_full_lattice_scale() {
    let tau = UpperHalfPlanePoint::tau_i();
    let truncation = LatticeSumTruncation::new(18).unwrap();
    let curve = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
    let report = curve
        .recover_period_basis(PeriodRecoveryConfig::strict())
        .unwrap();
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
    let report = curve
        .recover_canonical_tau(PeriodRecoveryConfig::strict())
        .unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        report.canonical_tau().tau(),
        tau.tau(),
        ApproxTolerance::new(1.0e-3, 1.0e-3)
    ));
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
            curve
                .recover_canonical_tau(PeriodRecoveryConfig::strict())
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
        curve in arb_stable_real_split_analytic_curve(),
    ) {
        let config = PeriodRecoveryConfig::strict();
        let report = curve.recover_canonical_tau(config)
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
