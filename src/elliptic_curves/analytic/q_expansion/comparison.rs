use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    ApproxTolerance, LatticeSumTruncation, QExpansionTruncation, UpperHalfPlanePoint,
};
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};

/// Structured comparison between the Eisenstein-sum and `q`-expansion routes
/// to the analytic modular `j`-invariant.
#[derive(Clone, Debug, PartialEq)]
pub struct JInvariantComparisonReport {
    tau: UpperHalfPlanePoint,
    comparison: ComplexApproxComparison,
    lattice_truncation: LatticeSumTruncation,
    q_truncation: QExpansionTruncation,
}

impl JInvariantComparisonReport {
    pub(crate) fn new(
        tau: UpperHalfPlanePoint,
        comparison: ComplexApproxComparison,
        lattice_truncation: LatticeSumTruncation,
        q_truncation: QExpansionTruncation,
    ) -> Self {
        Self {
            tau,
            comparison,
            lattice_truncation,
            q_truncation,
        }
    }

    /// Returns the upper-half-plane parameter `τ` used in both computations.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the `j`-invariant obtained from truncated Eisenstein sums.
    pub fn eisenstein_j(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the `j`-invariant obtained from the truncated `q`-expansion.
    pub fn q_expansion_j(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the residual `j_Eisenstein - j_q`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm `|j_Eisenstein - j_q|`.
    pub fn absolute_difference(&self) -> f64 {
        self.comparison.absolute_difference()
    }

    /// Returns whether the two approximations agreed under the supplied
    /// tolerance policy.
    pub fn agrees_approximately(&self) -> bool {
        self.comparison.agrees_approximately()
    }

    /// Returns the lattice-sum truncation used on the Eisenstein side.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.lattice_truncation
    }

    /// Returns the `q`-expansion truncation used on the cusp-expansion side.
    pub fn q_truncation(&self) -> QExpansionTruncation {
        self.q_truncation
    }

    /// Returns the comparison tolerance.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.comparison.tolerance()
    }
}

impl HasComplexApproxComparison for JInvariantComparisonReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}
