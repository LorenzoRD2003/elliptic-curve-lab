//! Approximate period-lattice metadata for future analytic period recovery.
//!
//! For a non-singular analytic Weierstrass model
//! `y² = 4x³ - g₂x - g₃`, one expects a complex lattice
//! `Λ = ℤω₁ + ℤω₂` whose associated torus `ℂ / Λ` uniformizes the curve.
//! The ordered basis is not unique: replacing `(ω₁, ω₂)` by another
//! positively oriented `SL₂(ℤ)`-equivalent basis describes the same lattice.
//!
//! This module therefore starts with small metadata objects rather than a
//! premature recovery algorithm. The current surface is meant to support:
//!
//! - storing one chosen approximate period basis
//! - recording the corresponding modulus `τ = ω₂ / ω₁`
//! - comparing the `j`-invariant implied by that recovered lattice against
//!   the original curve-side `j`
//!
//! TODO(period-recovery milestone):
//! implement an actual numerical recovery routine
//! `AnalyticWeierstrassCurve -> PeriodLatticeApprox`.
//! A mathematically natural path is:
//!
//! 1. solve for the cubic roots `e₁, e₂, e₃` of `4x³ - g₂x - g₃`
//! 2. recover a period basis from complete elliptic integrals or an AGM-based
//!    algorithm
//! 3. normalize the resulting basis to a preferred representative
//!    (for example by reducing `τ` to the standard fundamental domain)
//! 4. validate the recovery by comparing the curve-side and lattice-side
//!    `j`-invariants

use num_complex::Complex64;

use super::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ComplexApproxComparison, ComplexLattice,
    HasAnalyticLatticeContext, HasComplexApproxComparison, UpperHalfPlanePoint,
};
use crate::numerics::ApproxTolerance;

/// One chosen approximate period basis for the analytic uniformization
/// lattice of a complex elliptic curve.
///
/// The same curve admits many ordered bases related by `SL₂(ℤ)`, so this type
/// intentionally stores one chosen basis rather than pretending to represent a
/// canonical pair of periods. The cached modulus `τ = ω₂ / ω₁` is included
/// because it is the natural parameter for later modular-normalization and
/// `j`-comparison experiments.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodLatticeApprox {
    lattice: ComplexLattice,
    tau: UpperHalfPlanePoint,
}

impl PeriodLatticeApprox {
    /// Builds a period-lattice approximation from a validated lattice basis.
    pub fn new(lattice: ComplexLattice) -> Result<Self, AnalyticCurveError> {
        let tau = lattice.tau()?;
        Ok(Self { lattice, tau })
    }

    /// Builds the standard normalized basis `ω₁ = 1`, `ω₂ = τ`.
    ///
    /// This does not recover periods from a curve. It only packages the
    /// canonical lattice representative `Λ_τ = ℤ + ℤτ` when the modulus `τ`
    /// is already known.
    pub fn standard_from_tau(tau: UpperHalfPlanePoint) -> Self {
        Self {
            lattice: ComplexLattice::from_tau(tau.clone()),
            tau,
        }
    }

    /// Returns the stored lattice basis.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the first period `ω₁`.
    pub fn omega1(&self) -> &Complex64 {
        self.lattice.omega1()
    }

    /// Returns the second period `ω₂`.
    pub fn omega2(&self) -> &Complex64 {
        self.lattice.omega2()
    }

    /// Returns the associated modulus `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }
}

impl HasAnalyticLatticeContext for PeriodLatticeApprox {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}

/// Comparison report between a curve-side `j`-invariant and the `j` implied by
/// a recovered approximate period lattice.
///
/// This report does not prescribe how the lattice was recovered. It only
/// stores the recovered basis together with a residual-style comparison
/// `j_recovered - j_curve`.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodRecoveryReport {
    curve: AnalyticWeierstrassCurve,
    periods: PeriodLatticeApprox,
    comparison: ComplexApproxComparison,
}

impl PeriodRecoveryReport {
    /// Builds the report from a curve, one recovered period lattice, and the
    /// `j`-invariant computed on the recovery side.
    ///
    /// The caller supplies `recovered_j` because different future algorithms
    /// may produce it through different routes, for example from a recovered
    /// modulus `τ`, from lattice Eisenstein sums, or from another normalization
    /// procedure.
    pub fn new(
        curve: AnalyticWeierstrassCurve,
        periods: PeriodLatticeApprox,
        recovered_j: Complex64,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        let curve_j = curve.j_invariant()?;
        Ok(Self {
            curve,
            periods,
            comparison: ComplexApproxComparison::new(recovered_j, curve_j, tolerance),
        })
    }

    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the recovered approximate period lattice.
    pub fn periods(&self) -> &PeriodLatticeApprox {
        &self.periods
    }

    /// Returns the `j`-invariant produced on the recovery side.
    pub fn recovered_j(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the `j`-invariant computed directly from the curve.
    pub fn curve_j(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the residual `j_recovered - j_curve`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm of the residual.
    pub fn absolute_difference(&self) -> f64 {
        self.comparison.absolute_difference()
    }

    /// Returns whether the two `j`-invariants agreed approximately.
    pub fn agrees_approximately(&self) -> bool {
        self.comparison.agrees_approximately()
    }

    /// Returns the tolerance policy used by the comparison.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.comparison.tolerance()
    }
}

impl HasComplexApproxComparison for PeriodRecoveryReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}

impl HasAnalyticLatticeContext for PeriodRecoveryReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        self.periods.tau()
    }

    fn lattice(&self) -> &ComplexLattice {
        self.periods.lattice()
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{PeriodLatticeApprox, PeriodRecoveryReport};
    use crate::elliptic_curves::analytic::{
        AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice, HasAnalyticLatticeContext,
        HasComplexApproxComparison, LatticeSumTruncation, UpperHalfPlanePoint,
    };

    #[test]
    fn standard_from_tau_uses_the_standard_z_plus_z_tau_basis() {
        let tau = UpperHalfPlanePoint::tau_i();
        let periods = PeriodLatticeApprox::standard_from_tau(tau.clone());

        assert_eq!(periods.omega1(), &Complex64::new(1.0, 0.0));
        assert_eq!(periods.omega2(), tau.tau());
        assert_eq!(periods.tau(), &tau);
    }

    #[test]
    fn new_recovers_tau_from_the_supplied_lattice() {
        let lattice =
            ComplexLattice::new(Complex64::new(2.0, 0.0), Complex64::new(1.0, 2.0)).unwrap();
        let periods = PeriodLatticeApprox::new(lattice.clone()).unwrap();

        assert_eq!(periods.lattice(), &lattice);
        assert_eq!(periods.tau().tau(), &Complex64::new(0.5, 1.0));
    }

    #[test]
    fn recovery_report_compares_recovered_and_curve_side_j_values() {
        let tau = UpperHalfPlanePoint::tau_i();
        let curve =
            AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap())
                .unwrap();
        let periods = PeriodLatticeApprox::standard_from_tau(tau);
        let recovered_j = curve.j_invariant().unwrap();

        let report =
            PeriodRecoveryReport::new(curve, periods, recovered_j, ApproxTolerance::strict())
                .unwrap();

        assert_eq!(report.recovered_j(), report.curve_j());
        assert_eq!(report.difference(), &Complex64::new(0.0, 0.0));
        assert!(report.agrees_approximately());
    }

    #[test]
    fn recovery_report_reuses_the_shared_context_traits() {
        let tau = UpperHalfPlanePoint::tau_i();
        let curve =
            AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap())
                .unwrap();
        let periods = PeriodLatticeApprox::standard_from_tau(tau.clone());
        let report = PeriodRecoveryReport::new(
            curve,
            periods.clone(),
            Complex64::new(1728.0, 0.0),
            ApproxTolerance::loose(),
        )
        .unwrap();

        assert_eq!(report.tau(), periods.tau());
        assert_eq!(report.lattice(), periods.lattice());
        assert_eq!(report.left(), report.recovered_j());
        assert_eq!(report.right(), report.curve_j());
    }
}
