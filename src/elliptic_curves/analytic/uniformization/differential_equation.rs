use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, ComplexLattice, EllipticFunctionTruncation,
    LatticeSumTruncation, TorusToCurveValues,
};
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};

/// Outcome of checking the differential equation
/// `℘′(z)^2 = 4℘(z)^3 - g₂℘(z) - g₃`.
#[derive(Clone, Copy, Debug, PartialEq)]
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
    pub fn status(&self) -> WeierstrassDifferentialEquationStatus {
        self.status
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
impl ComplexLattice {
    /// Verifies the classical differential equation relating `℘`, `℘′`, and the
    /// analytic invariants attached to one lattice.
    ///
    /// Complexity: `Θ(r_inv² + r_fun²)`, where `r_inv` is the invariant
    /// truncation radius and `r_fun` is the elliptic-function truncation radius.
    pub fn verify_weierstrass_differential_equation(
        &self,
        z: Complex64,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<WeierstrassDifferentialEquationReport, AnalyticCurveError> {
        let map =
            self.map_torus_point_to_curve(z, invariant_truncation, function_truncation, tolerance)?;
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
}
