use num_complex::Complex64;

use super::super::{
    AnalyticCurveError, ApproxTolerance, LatticeSumTruncation, UpperHalfPlanePoint,
    analytic_invariants_from_tau,
};
use super::ModularMatrix;
use crate::fields::ComplexApprox;

/// Structured comparison between the truncated analytic values
/// `j(τ)` and `j(γτ)` for one modular transformation `γ`.
///
/// Mathematically the classical modular invariant satisfies `j(γτ) = j(τ)`
/// for every `γ ∈ SL_2(ℤ)`.
/// The reason is not that `γ` leaves the number `τ` unchanged, but that `τ`
/// and `γτ` describe the same complex torus up to isomorphism. In that sense,
/// `j` depends on the geometric torus represented by `τ`, not on the
/// particular upper-half-plane coordinate used to write it. Said differently,
/// `j` descends from `\mathfrak H` to the modular quotient
/// `SL_2(ℤ)\backslash\mathfrak H`.
///
/// This report is only a finite-truncation numerical experiment. Because the
/// current lattice sums use square boxes in basis coordinates, the finite
/// truncation itself is not modularly invariant term-by-term. So a visible
/// residual at fixed radius does not mean that `j` stopped being modular; it
/// can also mean that the two coordinate-dependent truncations have not yet
/// converged closely enough.
#[derive(Clone, Debug, PartialEq)]
pub struct ModularInvarianceReport {
    original_tau: UpperHalfPlanePoint,
    transformed_tau: UpperHalfPlanePoint,
    matrix: ModularMatrix,
    original_j: Complex64,
    transformed_j: Complex64,
    difference: Complex64,
    invariant_approximately: bool,
    truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
}

impl ModularInvarianceReport {
    /// Returns the original upper-half-plane input `τ`.
    pub fn original_tau(&self) -> &UpperHalfPlanePoint {
        &self.original_tau
    }

    /// Returns the transformed point `γτ`.
    pub fn transformed_tau(&self) -> &UpperHalfPlanePoint {
        &self.transformed_tau
    }

    /// Returns the modular matrix `γ`.
    pub fn matrix(&self) -> ModularMatrix {
        self.matrix
    }

    /// Returns the truncated analytic value `j(τ)`.
    pub fn original_j(&self) -> &Complex64 {
        &self.original_j
    }

    /// Returns the truncated analytic value `j(γτ)`.
    pub fn transformed_j(&self) -> &Complex64 {
        &self.transformed_j
    }

    /// Returns the residual `j(τ) - j(γτ)`.
    pub fn difference(&self) -> &Complex64 {
        &self.difference
    }

    /// Returns the Euclidean norm `|j(τ) - j(γτ)|`.
    pub fn absolute_difference(&self) -> f64 {
        self.difference.norm()
    }

    /// Returns whether the two truncated values agreed under the supplied
    /// tolerance policy.
    pub fn invariant_approximately(&self) -> bool {
        self.invariant_approximately
    }

    /// Returns the lattice-sum truncation used on both sides.
    pub fn truncation(&self) -> LatticeSumTruncation {
        self.truncation
    }

    /// Returns the tolerance used for the comparison verdict.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }
}

/// Compares the two truncated analytic values `j(τ)` and `j(γτ)`.
///
/// This first applies the modular transformation `γ` to `τ`, then computes the
/// analytic `j`-invariant on both upper-half-plane points using the same
/// lattice-sum truncation radius.
///
/// Conceptually this is checking whether our current numerical approximations
/// respect the fact that `τ` and `γτ` should encode the same underlying torus.
/// It is therefore a test of the approximation scheme as much as of the
/// mathematical identity itself.
///
/// Complexity: `Θ(r²)` in the truncation radius `r`, up to the constant factor
/// coming from evaluating the truncated lattice invariants at two points.
pub fn verify_j_modular_invariance(
    tau: UpperHalfPlanePoint,
    matrix: ModularMatrix,
    truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<ModularInvarianceReport, AnalyticCurveError> {
    let transformed_tau = matrix.apply(&tau)?;
    let original_invariants = analytic_invariants_from_tau(&tau, truncation)?;
    let transformed_invariants = analytic_invariants_from_tau(&transformed_tau, truncation)?;
    let comparison = ComplexApprox::comparison_report(
        &original_invariants.j_invariant,
        &transformed_invariants.j_invariant,
        tolerance,
    );

    Ok(ModularInvarianceReport {
        original_tau: tau,
        transformed_tau,
        matrix,
        original_j: original_invariants.j_invariant,
        transformed_j: transformed_invariants.j_invariant,
        difference: comparison.difference,
        invariant_approximately: comparison.is_close,
        truncation,
        tolerance,
    })
}
