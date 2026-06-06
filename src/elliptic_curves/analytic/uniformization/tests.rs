mod tests {
    use crate::elliptic_curves::analytic::{
        AnalyticCurvePoint, ApproxTolerance, ComplexLattice, EllipticFunctionApproximation,
        EllipticFunctionTruncation, FundamentalParallelogramCoordinate, LatticeSumTruncation,
        TorusToCurveValues, UpperHalfPlanePoint, WeierstrassDifferentialEquationStatus,
        map_fundamental_point_to_curve, map_torus_point_to_curve,
        verify_weierstrass_differential_equation, weierstrass_p, weierstrass_p_derivative,
    };
    use crate::fields::ComplexApprox;
    use num_complex::Complex64;
    use proptest::prelude::*;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn square_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
    }

    #[test]
    fn lattice_points_map_to_infinity_and_report_a_pole() {
        let lattice = square_lattice();
        let result = map_torus_point_to_curve(
            &lattice,
            c(1.0, 0.0),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::strict(),
        )
        .unwrap();

        assert_eq!(*result.z(), c(1.0, 0.0));
        assert_eq!(result.values(), &TorusToCurveValues::Pole);
        assert_eq!(result.point(), &AnalyticCurvePoint::infinity());
        assert!(result.lies_on_curve());
    }

    #[test]
    fn finite_torus_points_store_the_same_p_and_p_prime_as_direct_evaluation() {
        let lattice = square_lattice();
        let z = c(0.3, 0.2);
        let invariant_truncation = LatticeSumTruncation::larger_for_comparison();
        let function_truncation = EllipticFunctionTruncation::default_educational();
        let tolerance = ApproxTolerance::loose();

        let result = map_torus_point_to_curve(
            &lattice,
            z,
            invariant_truncation,
            function_truncation,
            tolerance,
        )
        .unwrap();
        let p = weierstrass_p(&lattice, z, function_truncation).unwrap();
        let p_prime = weierstrass_p_derivative(&lattice, z, function_truncation).unwrap();

        match result.values() {
            TorusToCurveValues::Pole => panic!("expected finite torus-to-curve values"),
            TorusToCurveValues::FiniteValues {
                p: stored_p,
                p_prime: stored_p_prime,
            } => {
                assert!(ComplexApprox::eq_with_tolerance(
                    stored_p,
                    p.value(),
                    tolerance
                ));
                assert!(ComplexApprox::eq_with_tolerance(
                    stored_p_prime,
                    p_prime.value(),
                    tolerance
                ));
            }
        }

        assert_eq!(result.membership_report().point(), result.point());
        assert!(result.membership_report().absolute_error().is_finite());
    }

    #[test]
    fn torus_point_maps_to_point_on_curve() {
        let lattice = square_lattice();
        let result = map_torus_point_to_curve(
            &lattice,
            c(0.3, 0.2),
            LatticeSumTruncation::new(12).unwrap(),
            EllipticFunctionTruncation::new(14).unwrap(),
            ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        assert!(matches!(
            result.values(),
            TorusToCurveValues::FiniteValues { .. }
        ));
        assert!(result.lies_on_curve());
        assert!(result.membership_report().absolute_error().is_finite());
    }

    #[test]
    fn mapping_a_fundamental_coordinate_matches_mapping_its_complex_representative() {
        let lattice = square_lattice();
        let coord = FundamentalParallelogramCoordinate::new(0.25, 0.4).unwrap();
        let z = lattice.point_from_fundamental_coordinates(coord.clone());
        let invariant_truncation = LatticeSumTruncation::default_educational();
        let function_truncation = EllipticFunctionTruncation::default_educational();
        let tolerance = ApproxTolerance::loose();

        let from_coord = map_fundamental_point_to_curve(
            &lattice,
            coord,
            invariant_truncation,
            function_truncation,
            tolerance,
        )
        .unwrap();
        let from_z = map_torus_point_to_curve(
            &lattice,
            z,
            invariant_truncation,
            function_truncation,
            tolerance,
        )
        .unwrap();

        assert_eq!(from_coord, from_z);
    }

    #[test]
    fn lattice_points_report_the_pole_case() {
        let lattice = square_lattice();
        let tolerance = ApproxTolerance::strict();

        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.0, 1.0),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            tolerance,
        )
        .unwrap();

        assert_eq!(report.values(), &TorusToCurveValues::Pole);
        assert_eq!(
            report.status(),
            &WeierstrassDifferentialEquationStatus::Pole
        );
        assert_eq!(report.lhs(), &c(0.0, 0.0));
        assert_eq!(report.rhs(), &c(0.0, 0.0));
        assert_eq!(report.difference(), &c(0.0, 0.0));
        assert!(!report.holds_approximately());
        assert_eq!(report.tolerance(), tolerance);
    }

    #[test]
    fn finite_report_exposes_same_tolerance_and_finite_residuals() {
        let lattice = square_lattice();
        let tolerance = ApproxTolerance::loose();

        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.3, 0.2),
            LatticeSumTruncation::larger_for_comparison(),
            EllipticFunctionTruncation::default_educational(),
            tolerance,
        )
        .unwrap();

        match report.values() {
            TorusToCurveValues::Pole => panic!("expected finite differential-equation values"),
            TorusToCurveValues::FiniteValues { .. } => {}
        }

        assert_eq!(*report.z(), c(0.3, 0.2));
        assert_eq!(report.tolerance(), tolerance);
        assert!(report.lhs().re.is_finite());
        assert!(report.lhs().im.is_finite());
        assert!(report.rhs().re.is_finite());
        assert!(report.rhs().im.is_finite());
        assert!(report.difference().re.is_finite());
        assert!(report.difference().im.is_finite());
    }

    #[test]
    fn finite_status_matches_the_residual_verdict() {
        let lattice = square_lattice();
        let strict = ApproxTolerance::strict();
        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.3, 0.2),
            LatticeSumTruncation::larger_for_comparison(),
            EllipticFunctionTruncation::default_educational(),
            strict,
        )
        .unwrap();

        let residual_is_small =
            ComplexApprox::eq_with_tolerance(report.lhs(), report.rhs(), strict);
        assert_eq!(report.holds_approximately(), residual_is_small);
        assert_eq!(
            report.status(),
            if residual_is_small {
                &WeierstrassDifferentialEquationStatus::HoldsApproximately
            } else {
                &WeierstrassDifferentialEquationStatus::FailsApproximately
            }
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn generic_fundamental_coordinates_match_their_complex_representatives(
            u in 0.0f64..1.0,
            v in 0.0f64..1.0,
        ) {
            let lattice = square_lattice();
            let coord = FundamentalParallelogramCoordinate::new(u, v).unwrap();
            let z = lattice.point_from_fundamental_coordinates(coord.clone());
            let invariant_truncation = LatticeSumTruncation::default_educational();
            let function_truncation = EllipticFunctionTruncation::default_educational();
            let tolerance = ApproxTolerance::loose();

            let from_coord = map_fundamental_point_to_curve(
                &lattice,
                coord,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();
            let from_z = map_torus_point_to_curve(
                &lattice,
                z,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();

            prop_assert_eq!(from_coord, from_z);
        }

        #[test]
        fn torus_to_curve_map_is_invariant_under_small_integer_lattice_shifts(
            u in 0.15f64..0.85,
            v in 0.15f64..0.85,
            m in -2i32..=2,
            n in -2i32..=2,
        ) {
            let lattice = square_lattice();
            let z = c(u, v);
            let shifted = z + c(m as f64, n as f64);
            let invariant_truncation = LatticeSumTruncation::larger_for_comparison();
            let function_truncation = EllipticFunctionTruncation::default_educational();
            let tolerance = ApproxTolerance::loose();
            let comparison_tolerance = ApproxTolerance::new(1.0e-9, 1.0e-9);

            let original = map_torus_point_to_curve(
                &lattice,
                z,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();
            let translated = map_torus_point_to_curve(
                &lattice,
                shifted,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();

            prop_assert_eq!(original.point().is_identity(), translated.point().is_identity());
            prop_assert_eq!(original.lies_on_curve(), translated.lies_on_curve());

            match (original.values(), translated.values()) {
                (
                    TorusToCurveValues::FiniteValues { p: p_left, p_prime: p_prime_left },
                    TorusToCurveValues::FiniteValues { p: p_right, p_prime: p_prime_right },
                ) => {
                    prop_assert!(ComplexApprox::eq_with_tolerance(
                        p_left,
                        p_right,
                        comparison_tolerance,
                    ));
                    prop_assert!(ComplexApprox::eq_with_tolerance(
                        p_prime_left,
                        p_prime_right,
                        comparison_tolerance,
                    ));
                    prop_assert!(ComplexApprox::eq_with_tolerance(
                        original.membership_report().lhs(),
                        translated.membership_report().lhs(),
                        comparison_tolerance,
                    ));
                    prop_assert!(ComplexApprox::eq_with_tolerance(
                        original.membership_report().rhs(),
                        translated.membership_report().rhs(),
                        comparison_tolerance,
                    ));
                    prop_assert!(ComplexApprox::eq_with_tolerance(
                        original.membership_report().difference(),
                        translated.membership_report().difference(),
                        comparison_tolerance,
                    ));
                }
                (TorusToCurveValues::Pole, TorusToCurveValues::Pole) => {}
                other => prop_assert!(false, "mismatched torus-to-curve cases: {other:?}"),
            }
        }

        #[test]
        fn generic_finite_differential_reports_match_the_residual_verdict(
            u in 0.15f64..0.85,
            v in 0.15f64..0.85,
        ) {
            let lattice = square_lattice();
            let tolerance = ApproxTolerance::strict();
            let report = verify_weierstrass_differential_equation(
                &lattice,
                c(u, v),
                LatticeSumTruncation::larger_for_comparison(),
                EllipticFunctionTruncation::default_educational(),
                tolerance,
            ).unwrap();

            match report.values() {
                TorusToCurveValues::Pole => panic!("sampled point should stay away from the lattice"),
                TorusToCurveValues::FiniteValues { .. } => {}
            }

            let residual_is_small = ComplexApprox::eq_with_tolerance(
                report.lhs(),
                report.rhs(),
                tolerance,
            );
            prop_assert_eq!(report.holds_approximately(), residual_is_small);
            prop_assert_eq!(
                report.status(),
                if residual_is_small {
                    &WeierstrassDifferentialEquationStatus::HoldsApproximately
                } else {
                    &WeierstrassDifferentialEquationStatus::FailsApproximately
                }
            );
        }
    }
}
