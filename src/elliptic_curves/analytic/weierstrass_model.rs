use core::fmt;

use num_complex::Complex64;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::fields::ComplexApprox;
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};
use crate::visualization::fields::format_complex;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, ComplexLattice, LatticeSumTruncation,
    UpperHalfPlanePoint, analytic_discriminant, analytic_invariants, analytic_invariants_from_tau,
    analytic_j_invariant,
};

/// Analytic curve points for the model `y² = 4x³ - g₂x - g₃`.
///
/// This is an alias to the crate's standard affine-point representation over
/// the approximate complex backend. Geometrically, the point shape is the same
/// as in the algebraic short-Weierstrass modules: either the distinguished
/// point at infinity or an affine pair `(x, y)`.
pub type AnalyticCurvePoint = AffinePoint<ComplexApprox>;

/// Structured approximate membership check for `y² = 4x³ - g₂x - g₃`.
///
/// For finite points, the report compares the left-hand side `y²` against the
/// right-hand side `4x³ - g₂x - g₃`. For the point at infinity, both sides are
/// recorded as zero and the report is marked as on-curve.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticCurveMembershipReport {
    /// Point being checked against the analytic Weierstrass equation.
    point: AnalyticCurvePoint,
    /// Shared complex approximate comparison payload for `y²` versus
    /// `4x³ - g₂x - g₃`.
    comparison: ComplexApproxComparison,
}

impl AnalyticCurveMembershipReport {
    /// Returns the checked point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        &self.point
    }

    /// Returns the left-hand side value, usually `y²`.
    pub fn lhs(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the right-hand side value, usually `4x³ - g₂x - g₃`.
    pub fn rhs(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the residual `lhs - rhs`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm `|lhs - rhs|`.
    pub fn absolute_error(&self) -> f64 {
        self.comparison.absolute_difference()
    }

    /// Returns the tolerance policy used by the comparison.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.comparison.tolerance()
    }

    /// Returns whether the point was accepted as lying on the curve.
    pub fn is_on_curve(&self) -> bool {
        self.comparison.agrees_approximately()
    }
}

impl HasComplexApproxComparison for AnalyticCurveMembershipReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}

/// Analytic Weierstrass model `y² = 4x³ - g₂x - g₃`
/// attached to approximate complex lattice invariants.
///
/// This wrapper keeps the analytic coefficients `g₂` and `g₃` as the public
/// surface, but internally it can be translated to the algebraic short
/// Weierstrass form `Y² = x³ + ax + b`, with `a = -g₂ / 4` and `b = -g₃ / 4`.
///
/// Under that change of variables `y = 2Y`, the discriminant condition is the
/// same, so non-singularity can reuse the existing short-Weierstrass
/// validation logic.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticWeierstrassCurve {
    g2: Complex64,
    g3: Complex64,
}

impl AnalyticWeierstrassCurve {
    /// Builds a validated analytic Weierstrass model from `g₂` and `g₃`.
    ///
    /// The constructor rejects when the  discriminant `Δ = g₂^3 - 27 g₃^2`
    /// is too close to zero for the approximate complex backend.
    pub fn new(g2: Complex64, g3: Complex64) -> Result<Self, AnalyticCurveError> {
        let curve = Self { g2, g3 };
        curve.try_as_short_weierstrass().map(|_| curve)
    }

    /// Builds the analytic Weierstrass model coming from one complex lattice.
    ///
    /// This uses the truncated analytic lattice invariants of `Λ` and then
    /// interprets them in the analytic equation `y² = 4x³ - g₂x - g₃`.
    ///
    /// Complexity: `Θ(r²)` in the invariant truncation radius `r`.
    pub fn from_lattice(
        lattice: &ComplexLattice,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let invariants = analytic_invariants(lattice, truncation)?;
        Self::new(invariants.g2, invariants.g3)
    }

    /// Builds the analytic Weierstrass model attached to the standard lattice
    /// `Λ_τ = ℤ + ℤτ`.
    ///
    /// Complexity: `Θ(r²)` in the invariant truncation radius `r`.
    pub fn from_tau(
        tau: &UpperHalfPlanePoint,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let invariants = analytic_invariants_from_tau(tau, truncation)?;
        Self::new(invariants.g2, invariants.g3)
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
        analytic_discriminant(&self.g2, &self.g3)
    }

    /// Returns the analytic `j`-invariant `j = 1728 g₂^3 / Δ`.
    pub fn j_invariant(&self) -> Result<Complex64, AnalyticCurveError> {
        analytic_j_invariant(&self.g2, &self.g3)
    }

    /// Returns the algebraic short-Weierstrass companion
    /// `Y² = x³ + ax + b` with `a = -g₂ / 4` and `b = -g₃ / 4`.
    ///
    /// This is the internal algebraic model corresponding to the analytic
    /// equation after the coordinate change `y = 2Y`.
    pub fn as_short_weierstrass(&self) -> ShortWeierstrassCurve<ComplexApprox> {
        self.try_as_short_weierstrass()
            .expect("validated analytic Weierstrass curve must convert to a non-singular short Weierstrass model")
    }

    /// Returns the right-hand side `4x³ - g₂x - g₃` of the analytic
    /// Weierstrass equation.
    pub fn rhs(&self, x: &Complex64) -> Complex64 {
        Complex64::new(4.0, 0.0) * x.powu(3) - self.g2 * *x - self.g3
    }

    /// Returns whether `point` approximately satisfies
    /// `y² = 4x³ - g₂x - g₃` under the caller-provided tolerance.
    ///
    /// The point at infinity is always accepted as belonging to the projective
    /// completion of the curve.
    pub fn contains_point_approx(
        &self,
        point: &AnalyticCurvePoint,
        tolerance: ApproxTolerance,
    ) -> bool {
        self.membership_report(point, tolerance).is_on_curve()
    }

    /// Returns a structured approximate membership report for `point`.
    ///
    /// For affine points, this compares `y²` against `4x³ - g₂x - g₃` under
    /// the supplied tolerance. For the point at infinity, the report is
    /// accepted automatically and records zero residual.
    pub fn membership_report(
        &self,
        point: &AnalyticCurvePoint,
        tolerance: ApproxTolerance,
    ) -> AnalyticCurveMembershipReport {
        match point {
            AffinePoint::Infinity => AnalyticCurveMembershipReport {
                point: point.clone(),
                comparison: ComplexApproxComparison::new(
                    Complex64::new(0.0, 0.0),
                    Complex64::new(0.0, 0.0),
                    tolerance,
                ),
            },
            AffinePoint::Finite { x, y } => {
                let left = y.powu(2);
                let right = self.rhs(x);
                AnalyticCurveMembershipReport {
                    point: point.clone(),
                    comparison: ComplexApproxComparison::new(left, right, tolerance),
                }
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

/// Small value object for the short-Weierstrass companion
/// `y² = x³ + ax + b` attached to an analytic cubic.
///
/// This keeps the short-form coefficients and their `j`-invariant visible
/// without exposing the full generic algebraic curve type inside a
/// higher-level educational report.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticShortWeierstrassModel {
    a: Complex64,
    b: Complex64,
    j_invariant: Complex64,
}

impl AnalyticShortWeierstrassModel {
    /// Builds the short-Weierstrass companion attached to an analytic cubic.
    pub fn from_analytic_curve(curve: &AnalyticWeierstrassCurve) -> Self {
        let short_curve: ShortWeierstrassCurve<ComplexApprox> = curve.as_short_weierstrass();

        Self {
            a: *short_curve.a(),
            b: *short_curve.b(),
            j_invariant: short_curve.j_invariant(),
        }
    }

    /// Returns the short-form coefficient `a`.
    pub fn a(&self) -> &Complex64 {
        &self.a
    }

    /// Returns the short-form coefficient `b`.
    pub fn b(&self) -> &Complex64 {
        &self.b
    }

    /// Returns the short-model `j`-invariant.
    pub fn j_invariant(&self) -> &Complex64 {
        &self.j_invariant
    }
}

impl fmt::Display for AnalyticShortWeierstrassModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "y^2 = x^3 + ({:?})x + ({:?})", self.a, self.b)
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::elliptic_curves::analytic::{AnalyticCurvePoint, AnalyticWeierstrassCurve};
    use crate::elliptic_curves::analytic::{
        AnalyticCurveError, ApproxTolerance, ComplexLattice, LatticeSumTruncation,
        UpperHalfPlanePoint, analytic_discriminant, analytic_invariants,
    };
    use crate::fields::{ComplexApprox, Field};

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn constructor_rejects_nearly_singular_coefficients() {
        assert_eq!(
            AnalyticWeierstrassCurve::new(c(0.0, 0.0), c(0.0, 0.0)),
            Err(AnalyticCurveError::NearlySingularAnalyticCurve)
        );
    }

    #[test]
    fn discriminant_matches_the_analytic_formula() {
        let curve = AnalyticWeierstrassCurve::new(c(12.0, 1.0), c(-3.0, 2.0)).unwrap();
        let expected = analytic_discriminant(curve.g2(), curve.g3());

        assert!(ComplexApprox::eq(&curve.discriminant(), &expected));
    }

    #[test]
    fn j_invariant_matches_the_shared_analytic_helper() {
        let curve = AnalyticWeierstrassCurve::new(c(12.0, 1.0), c(4.0, -2.0)).unwrap();
        let expected =
            crate::elliptic_curves::analytic::analytic_j_invariant(curve.g2(), curve.g3()).unwrap();

        assert!(ComplexApprox::eq(&curve.j_invariant().unwrap(), &expected));
    }

    #[test]
    fn rhs_uses_the_analytic_equation() {
        let curve = AnalyticWeierstrassCurve::new(c(8.0, 0.0), c(-4.0, 0.0)).unwrap();
        let x = c(2.0, 0.0);

        assert_eq!(
            curve.rhs(&x),
            c(4.0, 0.0) * x.powu(3) - c(8.0, 0.0) * x - c(-4.0, 0.0)
        );
    }

    #[test]
    fn approximate_membership_accepts_infinity_and_known_affine_points() {
        let curve = AnalyticWeierstrassCurve::new(c(0.0, 0.0), c(4.0, 0.0)).unwrap();
        let infinity = AnalyticCurvePoint::infinity();
        let affine = AnalyticCurvePoint::new(c(1.0, 0.0), c(0.0, 0.0));
        let off_curve = AnalyticCurvePoint::new(c(0.0, 0.0), c(0.0, 0.0));
        let tolerance = ApproxTolerance::strict();

        assert!(curve.contains_point_approx(&infinity, tolerance));
        assert!(curve.contains_point_approx(&affine, tolerance));
        assert!(!curve.contains_point_approx(&off_curve, tolerance));
    }

    #[test]
    fn membership_report_records_lhs_rhs_and_residual() {
        let curve = AnalyticWeierstrassCurve::new(c(0.0, 0.0), c(4.0, 0.0)).unwrap();
        let point = AnalyticCurvePoint::new(c(0.0, 0.0), c(0.0, 0.0));
        let tolerance = ApproxTolerance::strict();

        let report = curve.membership_report(&point, tolerance);

        assert_eq!(report.point(), &point);
        assert_eq!(report.lhs(), &c(0.0, 0.0));
        assert_eq!(report.rhs(), &c(-4.0, 0.0));
        assert_eq!(report.difference(), &c(4.0, 0.0));
        assert_eq!(report.absolute_error(), 4.0);
        assert_eq!(report.tolerance(), tolerance);
        assert!(!report.is_on_curve());
    }

    #[test]
    fn membership_report_accepts_infinity_with_zero_residual() {
        let curve = AnalyticWeierstrassCurve::new(c(3.0, 0.0), c(-2.0, 0.0)).unwrap();
        let tolerance = ApproxTolerance::strict();
        let report = curve.membership_report(&AnalyticCurvePoint::infinity(), tolerance);

        assert!(report.is_on_curve());
        assert_eq!(report.lhs(), &c(0.0, 0.0));
        assert_eq!(report.rhs(), &c(0.0, 0.0));
        assert_eq!(report.difference(), &c(0.0, 0.0));
        assert_eq!(report.absolute_error(), 0.0);
    }

    #[test]
    fn approximate_membership_can_accept_small_numerical_noise() {
        let curve = AnalyticWeierstrassCurve::new(c(0.0, 0.0), c(4.0, 0.0)).unwrap();
        let noisy = AnalyticCurvePoint::new(c(1.0, 0.0), c(1.0e-5, 0.0));

        assert!(curve.contains_point_approx(&noisy, ApproxTolerance::loose()));
        assert!(!curve.contains_point_approx(&noisy, ApproxTolerance::strict()));
    }

    #[test]
    fn equation_string_uses_the_analytic_equation_surface() {
        let curve = AnalyticWeierstrassCurve::new(c(3.0, -2.0), c(-1.0, 4.0)).unwrap();
        let equation = curve.equation_string();

        assert!(equation.contains("y^2 = 4x^3"));
        assert!(equation.contains("3.000000 - 2.000000i"));
        assert!(equation.contains("-1.000000 + 4.000000i"));
    }

    #[test]
    fn display_reuses_the_equation_surface() {
        let curve = AnalyticWeierstrassCurve::new(c(3.0, -2.0), c(-1.0, 4.0)).unwrap();

        assert_eq!(format!("{curve}"), curve.equation_string());
    }

    #[test]
    fn public_short_weierstrass_translation_uses_expected_coefficients() {
        let curve = AnalyticWeierstrassCurve::new(c(8.0, 4.0), c(-12.0, 20.0)).unwrap();
        let short = curve.as_short_weierstrass();

        assert!(ComplexApprox::eq(short.a(), &c(-2.0, -1.0)));
        assert!(ComplexApprox::eq(short.b(), &c(3.0, -5.0)));
    }

    #[test]
    fn from_lattice_uses_the_same_g2_and_g3_as_analytic_invariants() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let truncation = LatticeSumTruncation::default_educational();

        let curve = AnalyticWeierstrassCurve::from_lattice(&lattice, truncation).unwrap();
        let invariants = analytic_invariants(&lattice, truncation).unwrap();

        assert!(ComplexApprox::eq(curve.g2(), &invariants.g2));
        assert!(ComplexApprox::eq(curve.g3(), &invariants.g3));
    }

    #[test]
    fn from_tau_matches_from_lattice_for_standard_tau_lattice() {
        let tau = UpperHalfPlanePoint::tau_rho();
        let truncation = LatticeSumTruncation::larger_for_comparison();

        let from_tau = AnalyticWeierstrassCurve::from_tau(&tau, truncation).unwrap();
        let lattice = ComplexLattice::from_tau(tau);
        let from_lattice = AnalyticWeierstrassCurve::from_lattice(&lattice, truncation).unwrap();

        assert_eq!(from_tau, from_lattice);
    }

    #[test]
    fn short_model_j_matches_analytic_model_j() {
        let curve = AnalyticWeierstrassCurve::from_tau(
            &UpperHalfPlanePoint::tau_generic_example(),
            LatticeSumTruncation::new(12).unwrap(),
        )
        .unwrap();
        let short = curve.as_short_weierstrass();

        assert!(ComplexApprox::eq(
            &curve.j_invariant().unwrap(),
            &short.j_invariant()
        ));
    }
}
