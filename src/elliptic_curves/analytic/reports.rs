use num_complex::Complex64;

use super::{
    AnalyticCurveError, ApproxTolerance, ComplexLattice, EllipticFunctionTruncation,
    LatticeSumTruncation, TorusToCurveValues, map_torus_point_to_curve,
};

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
/// The finite case compares
/// `℘′(z)^2` against `4℘(z)^3 - g₂℘(z) - g₃` using the same curve-membership
/// report produced by the torus-to-curve map. In the pole case, the report
/// reuses the point-at-infinity convention from that membership report:
/// `lhs = rhs = difference = 0`, with status [`Pole`](WeierstrassDifferentialEquationStatus::Pole).
#[derive(Clone, Debug, PartialEq)]
pub struct WeierstrassDifferentialEquationReport {
    z: Complex64,
    values: TorusToCurveValues,
    lhs: Complex64,
    rhs: Complex64,
    difference: Complex64,
    status: WeierstrassDifferentialEquationStatus,
    tolerance: ApproxTolerance,
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
        &self.lhs
    }

    /// Returns the right-hand side of the verified equation.
    pub fn rhs(&self) -> &Complex64 {
        &self.rhs
    }

    /// Returns the residual `lhs - rhs`.
    pub fn difference(&self) -> &Complex64 {
        &self.difference
    }

    /// Returns the interpreted verification status.
    pub fn status(&self) -> &WeierstrassDifferentialEquationStatus {
        &self.status
    }

    /// Returns the tolerance used for the approximate comparison.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns whether the finite verification held approximately.
    pub fn holds_approximately(&self) -> bool {
        matches!(
            self.status,
            WeierstrassDifferentialEquationStatus::HoldsApproximately
        )
    }
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
        TorusToCurveValues::FiniteValues { .. } if membership.is_on_curve => {
            WeierstrassDifferentialEquationStatus::HoldsApproximately
        }
        TorusToCurveValues::FiniteValues { .. } => {
            WeierstrassDifferentialEquationStatus::FailsApproximately
        }
    };

    Ok(WeierstrassDifferentialEquationReport {
        z: *map.z(),
        values: map.values().clone(),
        lhs: membership.lhs,
        rhs: membership.rhs,
        difference: membership.difference,
        status,
        tolerance,
    })
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{WeierstrassDifferentialEquationStatus, verify_weierstrass_differential_equation};
    use crate::elliptic_curves::analytic::{
        ApproxTolerance, ComplexLattice, EllipticFunctionTruncation, LatticeSumTruncation,
        TorusToCurveValues, UpperHalfPlanePoint,
    };
    use crate::fields::ComplexApprox;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn square_lattice() -> ComplexLattice {
        ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
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
}
