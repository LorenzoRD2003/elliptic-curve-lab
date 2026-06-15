use num_complex::Complex64;

use crate::elliptic_curves::analytic::LatticeSumTruncation;
use crate::numerics::ComplexDifferenceReport;

/// Finite square-box approximation to the lattice Eisenstein sum `G_k(Λ)`.
#[derive(Clone, Debug, PartialEq)]
pub struct EisensteinSumApprox {
    weight: usize,
    value: Complex64,
    truncation: LatticeSumTruncation,
    terms_used: usize,
}

impl EisensteinSumApprox {
    pub(crate) fn new(
        weight: usize,
        value: Complex64,
        truncation: LatticeSumTruncation,
        terms_used: usize,
    ) -> Self {
        Self {
            weight,
            value,
            truncation,
            terms_used,
        }
    }

    /// Weight `k` in `G_k(Λ)`.
    pub fn weight(&self) -> usize {
        self.weight
    }

    /// Truncated lattice-sum value in `ℂ`.
    pub fn value(&self) -> &Complex64 {
        &self.value
    }

    /// Square-box truncation used for the approximation.
    pub fn truncation(&self) -> LatticeSumTruncation {
        self.truncation
    }

    /// Number of nonzero lattice terms conceptually included in the box.
    pub fn terms_used(&self) -> usize {
        self.terms_used
    }
}

/// Side-by-side comparison between two square-box truncations of the same
/// lattice Eisenstein sum.
#[derive(Clone, Debug, PartialEq)]
pub struct TruncationConvergenceReport {
    small: LatticeSumTruncation,
    large: LatticeSumTruncation,
    comparison: ComplexDifferenceReport,
}

impl TruncationConvergenceReport {
    pub(crate) fn new(
        small: LatticeSumTruncation,
        large: LatticeSumTruncation,
        comparison: ComplexDifferenceReport,
    ) -> Self {
        Self {
            small,
            large,
            comparison,
        }
    }

    /// Returns the smaller truncation radius.
    pub fn small(&self) -> LatticeSumTruncation {
        self.small
    }

    /// Returns the larger truncation radius.
    pub fn large(&self) -> LatticeSumTruncation {
        self.large
    }

    /// Returns the truncated value computed with `small`.
    pub fn small_value(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the truncated value computed with `large`.
    pub fn large_value(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the residual `large_value - small_value`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm `|large_value - small_value|`.
    pub fn absolute_difference(&self) -> f64 {
        self.comparison.absolute_difference()
    }
}
