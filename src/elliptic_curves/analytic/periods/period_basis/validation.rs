use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ComplexLattice, UpperHalfPlanePoint,
    lattice::HasAnalyticLatticeContext, periods::PeriodLatticeApprox,
};
use crate::numerics::{ApproxTolerance, ComplexApproxComparison, HasComplexApproxComparison};

/// Comparison report between a curve-side `j`-invariant and the `j` implied by
/// a recovered approximate period lattice.
///
/// This report belongs to the period-recovery layer rather than the
/// inverse-uniformization layer: it explains whether one chosen recovered
/// period lattice already lands in the correct modular class of the curve.
/// It does not prescribe how the lattice was recovered.
#[derive(Clone, Debug, PartialEq)]
pub struct CurvePeriodLatticeComparisonReport {
    curve: AnalyticWeierstrassCurve,
    periods: PeriodLatticeApprox,
    comparison: ComplexApproxComparison,
}

impl CurvePeriodLatticeComparisonReport {
    /// Builds the report from a curve, one recovered period lattice, and the
    /// `j`-invariant computed on the recovery side.
    ///
    /// The caller supplies `recovered_j` because different algorithms may
    /// obtain it through different routes, for example from a recovered
    /// modulus `τ`, from lattice Eisenstein sums, or from another
    /// normalization procedure.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn new(
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

impl HasComplexApproxComparison for CurvePeriodLatticeComparisonReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}

impl HasAnalyticLatticeContext for CurvePeriodLatticeComparisonReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        self.periods.tau()
    }

    fn lattice(&self) -> &ComplexLattice {
        self.periods.lattice()
    }
}
