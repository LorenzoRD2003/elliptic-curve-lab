use num_complex::Complex64;

use crate::elliptic_curves::analytic::{ApproxTolerance, weierstrass_model::AnalyticCurvePoint};
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};

/// Structured approximate membership check for `y² = 4x³ - g₂x - g₃`.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticCurveMembershipReport {
    point: AnalyticCurvePoint,
    comparison: ComplexApproxComparison,
}

impl AnalyticCurveMembershipReport {
    pub(crate) fn new(point: AnalyticCurvePoint, comparison: ComplexApproxComparison) -> Self {
        Self { point, comparison }
    }

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
