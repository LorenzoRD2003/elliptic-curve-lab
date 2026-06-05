use num_complex::Complex64;

use super::{
    AnalyticCurveError, AnalyticCurveMembershipReport, AnalyticCurvePoint,
    AnalyticWeierstrassCurve, ApproxTolerance, ComplexApproxComparison, ComplexLattice,
    EllipticFunctionApproximation, EllipticFunctionTruncation, FundamentalParallelogramCoordinate,
    HasComplexApproxComparison, LatticeSumTruncation, weierstrass_p, weierstrass_p_derivative,
};
use crate::fields::ComplexApprox;

/// Finite or pole-valued output of the torus-to-curve map `z ↦ (℘(z), ℘′(z))`.
///
/// At lattice points, the classical meromorphic coordinates `℘(z)` and
/// `℘′(z)` blow up. Geometrically, those torus points map to the point at
/// infinity on the corresponding analytic Weierstrass curve, so this enum
/// records either that pole case or the finite pair of evaluated coordinates.
#[derive(Clone, Debug, PartialEq)]
pub enum TorusToCurveValues {
    /// The input represents a lattice point, so the map lands at infinity.
    Pole,
    /// The input has finite `℘` and `℘′` values.
    FiniteValues { p: Complex64, p_prime: Complex64 },
}

/// Structured result of mapping one torus point to the analytic Weierstrass
/// curve attached to the same lattice.
#[derive(Clone, Debug, PartialEq)]
pub struct TorusToCurveMapResult {
    z: Complex64,
    point: AnalyticCurvePoint,
    curve: AnalyticWeierstrassCurve,
    values: TorusToCurveValues,
    membership_report: AnalyticCurveMembershipReport,
}

impl TorusToCurveMapResult {
    /// Returns the original complex torus representative requested by the caller.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns the resulting analytic curve point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        &self.point
    }

    /// Returns the analytic Weierstrass curve used for the map.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the finite `℘`/`℘′` pair or the pole case.
    pub fn values(&self) -> &TorusToCurveValues {
        &self.values
    }

    /// Returns the approximate curve-membership report for `point`.
    pub fn membership_report(&self) -> &AnalyticCurveMembershipReport {
        &self.membership_report
    }

    /// Returns whether the mapped point was accepted as lying on the curve
    /// under the caller-provided tolerance.
    pub fn lies_on_curve(&self) -> bool {
        self.membership_report.is_on_curve()
    }
}

/// Outcome of checking the differential equation
/// `℘′(z)^2 = 4℘(z)^3 - g₂℘(z) - g₃`.
#[derive(Clone, Debug, PartialEq)]
pub enum WeierstrassDifferentialEquationStatus {
    /// The finite approximations satisfied the differential equation within
    /// the requested tolerance.
    HoldsApproximately,
    /// The finite approximations did not satisfy the differential equation
    /// within the requested tolerance.
    FailsApproximately,
    /// The input represented a lattice point, so the torus-to-curve map lands
    /// at the point at infinity instead of finite `℘` and `℘′` values.
    Pole,
}

/// Structured verification report for the classical Weierstrass differential
/// equation.
///
/// The finite case compares `℘′(z)^2` against `4℘(z)^3 - g₂℘(z) - g₃` using
/// the same curve-membership report produced by the torus-to-curve map. In the
/// pole case, the report reuses the point-at-infinity convention from that
/// membership report: `lhs = rhs = difference = 0`, with status
/// [`Pole`](WeierstrassDifferentialEquationStatus::Pole).
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassDifferentialEquationReport {
    z: Complex64,
    values: TorusToCurveValues,
    comparison: ComplexApproxComparison,
    status: WeierstrassDifferentialEquationStatus,
}

impl WeierstrassDifferentialEquationReport {
    /// Returns the original evaluation point requested by the caller.
    pub fn z(&self) -> &Complex64 {
        &self.z
    }

    /// Returns the finite `℘`/`℘′` pair or the pole case.
    pub fn values(&self) -> &TorusToCurveValues {
        &self.values
    }

    /// Returns the left-hand side of the verified equation.
    pub fn lhs(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the right-hand side of the verified equation.
    pub fn rhs(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the residual `lhs - rhs`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the interpreted verification status.
    pub fn status(&self) -> &WeierstrassDifferentialEquationStatus {
        &self.status
    }

    /// Returns the tolerance used for the approximate comparison.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.comparison.tolerance()
    }

    /// Returns whether the finite verification held approximately.
    pub fn holds_approximately(&self) -> bool {
        matches!(
            self.status,
            WeierstrassDifferentialEquationStatus::HoldsApproximately
        )
    }
}

impl HasComplexApproxComparison for WeierstrassDifferentialEquationReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}

/// Maps one complex torus point to the analytic Weierstrass curve attached to
/// the same lattice.
///
/// This implements the analytic correspondence `z ↦ (℘(z), ℘′(z))`
/// between the complex torus `ℂ / Λ` and the cubic `y² = 4x³ - g₂x - g₃`.
///
/// If `z` represents a lattice point, then `℘(z)` and `℘′(z)` are poles and
/// the map lands at the point at infinity instead of returning an error.
///
/// Complexity: `Θ(r_inv² + r_fun²)`, where `r_inv` is the invariant
/// truncation radius and `r_fun` is the elliptic-function truncation radius.
pub fn map_torus_point_to_curve(
    lattice: &ComplexLattice,
    z: Complex64,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<TorusToCurveMapResult, AnalyticCurveError> {
    let reduced_coordinate = lattice.reduce_complex_point_to_fundamental_coordinates(z)?;
    let reduced_z = lattice.point_from_fundamental_coordinates(reduced_coordinate);
    let curve = AnalyticWeierstrassCurve::from_lattice(lattice, invariant_truncation)?;

    // reduced_z == 0 is equivalent to z being a lattice point
    if ComplexApprox::is_zero_with_tolerance(&reduced_z, ComplexApprox::default_tolerance()) {
        let point = AnalyticCurvePoint::infinity();
        let membership_report = curve.membership_report(&point, tolerance);

        return Ok(TorusToCurveMapResult {
            z,
            point,
            curve,
            values: TorusToCurveValues::Pole,
            membership_report,
        });
    }

    let p_approx = weierstrass_p(lattice, z, function_truncation)?;
    let p_prime_approx = weierstrass_p_derivative(lattice, z, function_truncation)?;
    let point = AnalyticCurvePoint::new(*p_approx.value(), *p_prime_approx.value());
    let membership_report = curve.membership_report(&point, tolerance);

    Ok(TorusToCurveMapResult {
        z,
        point: point.clone(),
        curve,
        values: TorusToCurveValues::FiniteValues {
            p: *p_approx.value(),
            p_prime: *p_prime_approx.value(),
        },
        membership_report,
    })
}

/// Maps one canonical fundamental-parallelogram coordinate to the analytic
/// Weierstrass curve attached to the same lattice.
///
/// Complexity: `Θ(r_inv² + r_fun²)`, inherited from [`map_torus_point_to_curve`].
pub fn map_fundamental_point_to_curve(
    lattice: &ComplexLattice,
    coord: FundamentalParallelogramCoordinate,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<TorusToCurveMapResult, AnalyticCurveError> {
    let z = lattice.point_from_fundamental_coordinates(coord);
    map_torus_point_to_curve(
        lattice,
        z,
        invariant_truncation,
        function_truncation,
        tolerance,
    )
}

/// Verifies the classical differential equation relating `℘`, `℘′`, and the
/// analytic invariants attached to one lattice.
///
/// Complexity: `Θ(r_inv² + r_fun²)`, where `r_inv` is the invariant
/// truncation radius and `r_fun` is the elliptic-function truncation radius.
pub fn verify_weierstrass_differential_equation(
    lattice: &ComplexLattice,
    z: Complex64,
    invariant_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
) -> Result<WeierstrassDifferentialEquationReport, AnalyticCurveError> {
    let map = map_torus_point_to_curve(
        lattice,
        z,
        invariant_truncation,
        function_truncation,
        tolerance,
    )?;
    let membership = map.membership_report();
    let status = match map.values() {
        TorusToCurveValues::Pole => WeierstrassDifferentialEquationStatus::Pole,
        TorusToCurveValues::FiniteValues { .. } if membership.is_on_curve() => {
            WeierstrassDifferentialEquationStatus::HoldsApproximately
        }
        TorusToCurveValues::FiniteValues { .. } => {
            WeierstrassDifferentialEquationStatus::FailsApproximately
        }
    };

    Ok(WeierstrassDifferentialEquationReport {
        z: *map.z(),
        values: map.values().clone(),
        comparison: match map.values() {
            TorusToCurveValues::Pole => ComplexApproxComparison::from_values_and_verdict(
                *membership.lhs(),
                *membership.rhs(),
                tolerance,
                false,
            ),
            TorusToCurveValues::FiniteValues { .. } => {
                ComplexApproxComparison::from_values_and_verdict(
                    *membership.lhs(),
                    *membership.rhs(),
                    tolerance,
                    membership.is_on_curve(),
                )
            }
        },
        status,
    })
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;
    use proptest::prelude::*;

    use super::{
        TorusToCurveValues, WeierstrassDifferentialEquationStatus, map_fundamental_point_to_curve,
        map_torus_point_to_curve, verify_weierstrass_differential_equation,
    };
    use crate::elliptic_curves::analytic::{
        AnalyticCurvePoint, ApproxTolerance, ComplexLattice, EllipticFunctionApproximation,
        EllipticFunctionTruncation, FundamentalParallelogramCoordinate, LatticeSumTruncation,
        UpperHalfPlanePoint, weierstrass_p, weierstrass_p_derivative,
    };
    use crate::fields::ComplexApprox;

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
