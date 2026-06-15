use num_complex::Complex64;
use proptest::prelude::*;

use crate::elliptic_curves::analytic::uniformization::UniformizationExperimentReport;
use crate::elliptic_curves::analytic::{
    AnalyticCurvePoint, ApproxTolerance, ComplexLattice, EllipticFunctionTruncation,
    FundamentalParallelogramCoordinate, LatticeSumTruncation, TorusToCurveValues,
    UpperHalfPlanePoint, WeierstrassDifferentialEquationStatus,
};
use crate::fields::complex_approx::ComplexApprox;
use crate::proptest_support::elliptic_curves::{
    arb_fundamental_coordinate, arb_interior_fundamental_coordinate,
};

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

#[test]
fn lattice_points_map_to_infinity_and_report_a_pole() {
    let lattice = square_lattice();
    let result = lattice
        .map_torus_point_to_curve(
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

    let result = lattice
        .map_torus_point_to_curve(z, invariant_truncation, function_truncation, tolerance)
        .unwrap();
    let p = lattice.weierstrass_p(z, function_truncation).unwrap();
    let p_prime = lattice
        .weierstrass_p_derivative(z, function_truncation)
        .unwrap();

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
    let result = lattice
        .map_torus_point_to_curve(
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

    let from_coord = lattice
        .map_fundamental_point_to_curve(coord, invariant_truncation, function_truncation, tolerance)
        .unwrap();
    let from_z = lattice
        .map_torus_point_to_curve(z, invariant_truncation, function_truncation, tolerance)
        .unwrap();

    assert_eq!(from_coord, from_z);
}

#[test]
fn lattice_points_report_the_pole_case() {
    let lattice = square_lattice();
    let tolerance = ApproxTolerance::strict();

    let report = lattice
        .verify_weierstrass_differential_equation(
            c(0.0, 1.0),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            tolerance,
        )
        .unwrap();

    assert_eq!(report.values(), &TorusToCurveValues::Pole);
    assert_eq!(report.status(), WeierstrassDifferentialEquationStatus::Pole);
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

    let report = lattice
        .verify_weierstrass_differential_equation(
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
    let report = lattice
        .verify_weierstrass_differential_equation(
            c(0.3, 0.2),
            LatticeSumTruncation::larger_for_comparison(),
            EllipticFunctionTruncation::default_educational(),
            strict,
        )
        .unwrap();

    let residual_is_small = ComplexApprox::eq_with_tolerance(report.lhs(), report.rhs(), strict);
    assert_eq!(report.holds_approximately(), residual_is_small);
    assert_eq!(
        report.status(),
        if residual_is_small {
            WeierstrassDifferentialEquationStatus::HoldsApproximately
        } else {
            WeierstrassDifferentialEquationStatus::FailsApproximately
        }
    );
}

#[test]
fn uniformization_report_derives_global_curve_membership_from_samples() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![
            Complex64::new(0.0, 0.0),
            Complex64::new(0.3, 0.2),
            Complex64::new(0.5, 0.0),
        ],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert_eq!(
        report.all_points_lie_on_curve(),
        report
            .sampled_points()
            .iter()
            .all(|point| point.lies_on_curve())
    );
}

#[test]
fn uniformization_report_keeps_the_same_curve_across_all_samples() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![Complex64::new(0.3, 0.2), Complex64::new(0.5, 0.0)],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert!(
        report
            .sampled_points()
            .iter()
            .all(|point| point.curve() == report.curve())
    );
}

#[test]
fn uniformization_report_can_include_both_finite_points_and_a_pole() {
    let report = UniformizationExperimentReport::from_sample_points(
        UpperHalfPlanePoint::tau_i(),
        vec![Complex64::new(0.0, 0.0), Complex64::new(0.3, 0.2)],
        LatticeSumTruncation::new(16).unwrap(),
        EllipticFunctionTruncation::new(14).unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2),
    )
    .unwrap();

    assert!(matches!(
        report.sampled_points()[0].values(),
        TorusToCurveValues::Pole
    ));
    assert!(matches!(
        report.sampled_points()[1].values(),
        TorusToCurveValues::FiniteValues { .. }
    ));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn generic_fundamental_coordinates_match_their_complex_representatives(
            coord in arb_fundamental_coordinate(),
        ) {
            let lattice = square_lattice();
            let z = lattice.point_from_fundamental_coordinates(coord.clone());
            let invariant_truncation = LatticeSumTruncation::default_educational();
            let function_truncation = EllipticFunctionTruncation::default_educational();
            let tolerance = ApproxTolerance::loose();

            let from_coord = lattice.map_fundamental_point_to_curve(
                coord,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();
            let from_z = lattice.map_torus_point_to_curve(
                z,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();

            prop_assert_eq!(from_coord, from_z);
        }

        #[test]
        fn torus_to_curve_map_is_invariant_under_small_integer_lattice_shifts(
            coord in arb_interior_fundamental_coordinate(),
            m in -2i32..=2,
            n in -2i32..=2,
        ) {
            let lattice = square_lattice();
            let z = lattice.point_from_fundamental_coordinates(coord);
            let shifted = z + c(m as f64, n as f64);
            let invariant_truncation = LatticeSumTruncation::larger_for_comparison();
            let function_truncation = EllipticFunctionTruncation::default_educational();
            let tolerance = ApproxTolerance::loose();
            let comparison_tolerance = ApproxTolerance::new(1.0e-9, 1.0e-9);

            let original = lattice.map_torus_point_to_curve(
                z,
                invariant_truncation,
                function_truncation,
                tolerance,
            ).unwrap();
            let translated = lattice.map_torus_point_to_curve(
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
            coord in arb_interior_fundamental_coordinate(),
        ) {
            let lattice = square_lattice();
            let tolerance = ApproxTolerance::strict();
            let report = lattice.verify_weierstrass_differential_equation(
                lattice.point_from_fundamental_coordinates(coord),
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
                    WeierstrassDifferentialEquationStatus::HoldsApproximately
                } else {
                    WeierstrassDifferentialEquationStatus::FailsApproximately
                }
            );
        }
}
