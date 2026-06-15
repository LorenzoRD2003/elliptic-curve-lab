use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
    LatticeSumTruncation, UpperHalfPlanePoint, lattice::HasAnalyticLatticeContext,
};
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};

/// Validation report for the inverse-uniformization direction
/// `τ -> Λ_τ -> (g₂, g₃, Δ, j)`.
///
/// Starting from one candidate upper-half-plane parameter `τ`, this report
/// recomputes the standard lattice `Λ_τ = ℤ + ℤτ`, approximates its analytic
/// invariants by a finite square-box lattice sum, and compares the recovered
/// modular `j`-invariant against the `j` attached directly to the target
/// analytic curve.
///
/// This is intentionally only a `j`-level validation. Agreement of `j`
/// suggests that the recovered `τ` describes a torus in the same modular
/// isomorphism class as the curve, but it does not by itself certify that all
/// scale-sensitive quantities such as `g₂` and `g₃` match in the same chosen
/// normalization.
#[derive(Clone, Debug, PartialEq)]
pub struct InverseUniformizationJValidationReport {
    curve: AnalyticWeierstrassCurve,
    tau: UpperHalfPlanePoint,
    lattice: ComplexLattice,
    recovered_invariants: AnalyticInvariants,
    comparison: ComplexApproxComparison,
}

impl InverseUniformizationJValidationReport {
    pub(crate) fn new(
        curve: AnalyticWeierstrassCurve,
        tau: UpperHalfPlanePoint,
        lattice: ComplexLattice,
        recovered_invariants: AnalyticInvariants,
        comparison: ComplexApproxComparison,
    ) -> Self {
        Self {
            curve,
            tau,
            lattice,
            recovered_invariants,
            comparison,
        }
    }

    /// Returns the target analytic curve whose `j`-invariant we validate
    /// against.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the candidate upper-half-plane parameter `τ`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the standard lattice `Λ_τ = ℤ + ℤτ`.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the analytic invariants recovered from `τ`.
    pub fn recovered_invariants(&self) -> &AnalyticInvariants {
        &self.recovered_invariants
    }

    /// Returns the recomputed lattice-side value `j(Λ_τ)`.
    pub fn recovered_j(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the curve-side value `j(E)`.
    pub fn curve_j(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the residual `j(Λ_τ) - j(E)`.
    pub fn difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm `|j(Λ_τ) - j(E)|`.
    pub fn absolute_difference(&self) -> f64 {
        self.comparison.absolute_difference()
    }

    /// Returns whether the two `j`-values agreed under the supplied tolerance.
    pub fn agrees_approximately(&self) -> bool {
        self.comparison.agrees_approximately()
    }

    /// Returns the lattice-sum truncation used to recover the analytic
    /// invariants from `τ`.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.recovered_invariants.truncation()
    }

    /// Returns the tolerance policy used for the residual verdict.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.comparison.tolerance()
    }
}

impl HasComplexApproxComparison for InverseUniformizationJValidationReport {
    fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }
}

impl HasAnalyticLatticeContext for InverseUniformizationJValidationReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}
