use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ComplexLattice, HasAnalyticLatticeContext,
    PeriodLatticeApprox, UpperHalfPlanePoint,
};
use crate::numerics::{ApproxTolerance, ComplexApproxComparison, HasComplexApproxComparison};

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
