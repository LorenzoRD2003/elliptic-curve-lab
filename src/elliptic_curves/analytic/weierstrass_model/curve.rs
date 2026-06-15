use core::fmt;

use num_complex::Complex64;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    analytic::{
        AnalyticCurveError, AnalyticInvariants, ApproxTolerance, ComplexLattice,
        LatticeSumTruncation, UpperHalfPlanePoint,
        weierstrass_model::membership::AnalyticCurveMembershipReport,
    },
};
use crate::fields::complex_approx::ComplexApprox;
use crate::numerics::ComplexApproxComparison;
use crate::visualization::fields::format_complex;

/// Analytic curve points for the model `y² = 4x³ - g₂x - g₃`.
pub type AnalyticCurvePoint = AffinePoint<ComplexApprox>;

/// Analytic Weierstrass model `y² = 4x³ - g₂x - g₃`
/// attached to approximate complex lattice invariants.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticWeierstrassCurve {
    g2: Complex64,
    g3: Complex64,
}

impl AnalyticWeierstrassCurve {
    /// Builds a validated analytic Weierstrass model from `g₂` and `g₃`.
    pub fn new(g2: Complex64, g3: Complex64) -> Result<Self, AnalyticCurveError> {
        let curve = Self { g2, g3 };
        curve.try_as_short_weierstrass().map(|_| curve)
    }

    /// Builds the analytic Weierstrass model coming from one complex lattice.
    pub fn from_lattice(
        lattice: &ComplexLattice,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let invariants = lattice.analytic_invariants(truncation)?;
        Self::new(*invariants.g2(), *invariants.g3())
    }

    /// Builds the analytic Weierstrass model attached to the standard lattice
    /// `Λ_τ = ℤ + ℤτ`.
    pub fn from_tau(
        tau: &UpperHalfPlanePoint,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let invariants = tau.analytic_invariants(truncation)?;
        Self::new(*invariants.g2(), *invariants.g3())
    }

    /// Returns the stored analytic invariant `g₂`.
    pub fn g2(&self) -> &Complex64 {
        &self.g2
    }

    /// Returns the stored analytic invariant `g₃`.
    pub fn g3(&self) -> &Complex64 {
        &self.g3
    }

    /// Returns the analytic discriminant `Δ = g₂^3 - 27 g₃^2`.
    pub fn discriminant(&self) -> Complex64 {
        AnalyticInvariants::discriminant_from_g2_g3(&self.g2, &self.g3)
    }

    /// Returns the analytic `j`-invariant `j = 1728 g₂^3 / Δ`.
    pub fn j_invariant(&self) -> Result<Complex64, AnalyticCurveError> {
        AnalyticInvariants::j_invariant_from_g2_g3(&self.g2, &self.g3)
    }

    /// Returns the algebraic short-Weierstrass companion.
    pub fn as_short_weierstrass(&self) -> ShortWeierstrassCurve<ComplexApprox> {
        self.try_as_short_weierstrass()
            .expect("validated analytic Weierstrass curve must convert to a non-singular short Weierstrass model")
    }

    /// Returns the right-hand side `4x³ - g₂x - g₃` of the analytic equation.
    pub fn rhs(&self, x: &Complex64) -> Complex64 {
        Complex64::new(4.0, 0.0) * x.powu(3) - self.g2 * *x - self.g3
    }

    /// Returns whether `point` approximately satisfies the analytic equation.
    pub fn contains_point_approx(
        &self,
        point: &AnalyticCurvePoint,
        tolerance: ApproxTolerance,
    ) -> bool {
        self.membership_report(point, tolerance).is_on_curve()
    }

    /// Returns a structured approximate membership report for `point`.
    pub fn membership_report(
        &self,
        point: &AnalyticCurvePoint,
        tolerance: ApproxTolerance,
    ) -> AnalyticCurveMembershipReport {
        match point {
            AffinePoint::Infinity => AnalyticCurveMembershipReport::new(
                point.clone(),
                ComplexApproxComparison::new(
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    tolerance,
                ),
            ),
            AffinePoint::Finite { x, y } => {
                let left = y.powu(2);
                let right = self.rhs(x);
                AnalyticCurveMembershipReport::new(
                    point.clone(),
                    ComplexApproxComparison::new(left, right, tolerance),
                )
            }
        }
    }

    /// Returns a human-readable equation string for the analytic model.
    pub fn equation_string(&self) -> String {
        format!(
            "y^2 = 4x^3 - ({})x - ({})",
            format_complex(&self.g2),
            format_complex(&self.g3)
        )
    }

    fn try_as_short_weierstrass(
        &self,
    ) -> Result<ShortWeierstrassCurve<ComplexApprox>, AnalyticCurveError> {
        ShortWeierstrassCurve::new(self.short_a(), self.short_b()).map_err(Self::map_curve_error)
    }

    fn short_a(&self) -> Complex64 {
        -self.g2 / Complex64::new(4.0, 0.0)
    }

    fn short_b(&self) -> Complex64 {
        -self.g3 / Complex64::new(4.0, 0.0)
    }

    fn map_curve_error(error: CurveError) -> AnalyticCurveError {
        match error {
            CurveError::SingularCurve => AnalyticCurveError::NearlySingularAnalyticCurve,
            _ => AnalyticCurveError::NumericalComparisonFailed,
        }
    }
}

impl fmt::Display for AnalyticWeierstrassCurve {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.equation_string())
    }
}
