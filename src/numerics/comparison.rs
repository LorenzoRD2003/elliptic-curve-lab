use num_complex::Complex64;

use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Exact left/right bookkeeping for a complex residual-style comparison.
///
/// This small value object keeps the two compared complex values together with
/// the residual `left - right`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ComplexDifferenceReport {
    left: Complex64,
    right: Complex64,
    difference: Complex64,
}

impl ComplexDifferenceReport {
    /// Builds the report from two compared complex values.
    pub fn new(left: Complex64, right: Complex64) -> Self {
        Self {
            left,
            right,
            difference: left - right,
        }
    }

    /// Returns the left-hand value.
    pub fn left(&self) -> &Complex64 {
        &self.left
    }

    /// Returns the right-hand value.
    pub fn right(&self) -> &Complex64 {
        &self.right
    }

    /// Returns the residual `left - right`.
    pub fn difference(&self) -> &Complex64 {
        &self.difference
    }

    /// Returns the Euclidean norm `|left - right|`.
    pub fn absolute_difference(&self) -> f64 {
        self.difference.norm()
    }
}

/// Approximate comparison between two complex values under one tolerance
/// policy.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ComplexApproxComparison {
    values: ComplexDifferenceReport,
    tolerance: ApproxTolerance,
    agrees_approximately: bool,
}

impl ComplexApproxComparison {
    /// Builds the comparison and derives the approximate verdict from the
    /// supplied tolerance.
    pub fn new(left: Complex64, right: Complex64, tolerance: ApproxTolerance) -> Self {
        Self::from_values_and_verdict(
            left,
            right,
            tolerance,
            ComplexApprox::eq_with_tolerance(&left, &right, tolerance),
        )
    }

    /// Builds the comparison while letting the caller override the approximate
    /// verdict.
    ///
    /// This is useful when the compared values are still worth recording but
    /// the mathematical meaning of the overall report requires a different
    /// status, for example a pole case that should not count as “holding
    /// approximately” even when both stored sides are zero by convention.
    pub fn from_values_and_verdict(
        left: Complex64,
        right: Complex64,
        tolerance: ApproxTolerance,
        agrees_approximately: bool,
    ) -> Self {
        Self {
            values: ComplexDifferenceReport::new(left, right),
            tolerance,
            agrees_approximately,
        }
    }

    /// Returns the compared values together with their residual.
    #[allow(dead_code)]
    pub fn values(&self) -> &ComplexDifferenceReport {
        &self.values
    }

    /// Returns the left-hand value.
    pub fn left(&self) -> &Complex64 {
        self.values.left()
    }

    /// Returns the right-hand value.
    pub fn right(&self) -> &Complex64 {
        self.values.right()
    }

    /// Returns the residual `left - right`.
    pub fn difference(&self) -> &Complex64 {
        self.values.difference()
    }

    /// Returns the Euclidean norm `|left - right|`.
    pub fn absolute_difference(&self) -> f64 {
        self.values.absolute_difference()
    }

    /// Returns the tolerance used for the approximate verdict.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns whether the two values were accepted as approximately equal.
    pub fn agrees_approximately(&self) -> bool {
        self.agrees_approximately
    }
}

/// Shared access to an embedded complex approximate comparison.
#[allow(dead_code)]
pub(crate) trait HasComplexApproxComparison {
    /// Returns the underlying comparison payload.
    fn comparison(&self) -> &ComplexApproxComparison;

    /// Returns the compared left-hand value.
    fn left(&self) -> &Complex64 {
        self.comparison().left()
    }

    /// Returns the compared right-hand value.
    fn right(&self) -> &Complex64 {
        self.comparison().right()
    }

    /// Returns the residual `left - right`.
    fn difference(&self) -> &Complex64 {
        self.comparison().difference()
    }

    /// Returns the Euclidean norm of the residual.
    fn absolute_difference(&self) -> f64 {
        self.comparison().absolute_difference()
    }

    /// Returns the tolerance used by the comparison.
    fn tolerance(&self) -> ApproxTolerance {
        self.comparison().tolerance()
    }

    /// Returns whether the compared values agreed approximately.
    fn agrees_approximately(&self) -> bool {
        self.comparison().agrees_approximately()
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::numerics::ApproxTolerance;
    use crate::numerics::{
        ComplexApproxComparison, ComplexDifferenceReport, HasComplexApproxComparison,
    };

    #[test]
    fn difference_report_stores_both_values_and_their_residual() {
        let report =
            ComplexDifferenceReport::new(Complex64::new(3.0, -1.0), Complex64::new(1.5, 2.0));

        assert_eq!(report.left(), &Complex64::new(3.0, -1.0));
        assert_eq!(report.right(), &Complex64::new(1.5, 2.0));
        assert_eq!(report.difference(), &Complex64::new(1.5, -3.0));
    }

    #[test]
    fn approximate_comparison_can_compute_its_own_verdict() {
        let comparison = ComplexApproxComparison::new(
            Complex64::new(1.0, 0.0),
            Complex64::new(1.0 + 1.0e-6, 0.0),
            ApproxTolerance::new(1.0e-3, 1.0e-3),
        );

        assert!(comparison.agrees_approximately());
        assert!((comparison.absolute_difference() - 1.0e-6).abs() <= 1.0e-12);
    }

    #[test]
    fn approximate_comparison_can_override_the_verdict() {
        let comparison = ComplexApproxComparison::from_values_and_verdict(
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            ApproxTolerance::strict(),
            false,
        );

        assert_eq!(comparison.difference(), &Complex64::new(0.0, 0.0));
        assert!(!comparison.agrees_approximately());
    }

    #[test]
    fn trait_exposes_shared_accessors() {
        #[derive(Clone, Debug, PartialEq)]
        struct Wrapper {
            comparison: ComplexApproxComparison,
        }

        impl HasComplexApproxComparison for Wrapper {
            fn comparison(&self) -> &ComplexApproxComparison {
                &self.comparison
            }
        }

        let wrapper = Wrapper {
            comparison: ComplexApproxComparison::new(
                Complex64::new(2.0, 0.0),
                Complex64::new(1.0, 0.0),
                ApproxTolerance::strict(),
            ),
        };

        assert_eq!(wrapper.left(), &Complex64::new(2.0, 0.0));
        assert_eq!(wrapper.right(), &Complex64::new(1.0, 0.0));
        assert_eq!(wrapper.difference(), &Complex64::new(1.0, 0.0));
        assert_eq!(wrapper.absolute_difference(), 1.0);
        assert_eq!(wrapper.tolerance(), ApproxTolerance::strict());
        assert!(!wrapper.agrees_approximately());
    }
}
