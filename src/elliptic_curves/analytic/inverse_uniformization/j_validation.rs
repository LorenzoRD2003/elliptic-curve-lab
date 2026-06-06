use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance,
    CanonicalTauRecoveryReport, ComplexLattice, HasAnalyticLatticeContext, LatticeSumTruncation,
    TauRecoveryReport, UpperHalfPlanePoint,
};
use crate::numerics::{ComplexApproxComparison, HasComplexApproxComparison};

use crate::elliptic_curves::analytic::inverse_uniformization::validation_shared::{
    compare_recovered_invariants_against_curve, recover_invariant_snapshot_from_tau,
};

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
        self.recovered_invariants.truncation
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

/// Validates one recovered upper-half-plane parameter `τ` against a target
/// analytic curve by comparing modular `j`-invariants.
///
/// The algorithm is:
///
/// 1. form the standard lattice `Λ_τ = ℤ + ℤτ`,
/// 2. approximate `g₂(Λ_τ), g₃(Λ_τ), Δ(Λ_τ), j(Λ_τ)` by finite lattice sums,
/// 3. compute the curve-side `j(E)` from the supplied Weierstrass model,
/// 4. compare `j(Λ_τ)` and `j(E)` under the requested tolerance.
///
/// The comparison is mathematically meaningful because `j` is the homothety-
/// and modular-invariant quantity shared by analytically uniformized and
/// algebraically presented elliptic curves. At the same time, the recovered
/// `j(Λ_τ)` is only approximate because the current `g₂` and `g₃` are
/// themselves computed from finite square-box truncations.
///
/// Complexity: `Θ(r²)` in the lattice truncation radius `r`, since the
/// dominant work is recomputing the truncated Eisenstein sums behind
/// `g₂(Λ_τ)` and `g₃(Λ_τ)`.
pub fn validate_recovered_tau_by_j_invariant(
    curve: &AnalyticWeierstrassCurve,
    tau: &UpperHalfPlanePoint,
    lattice_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
    let snapshot = recover_invariant_snapshot_from_tau(tau, lattice_truncation)?;
    let comparisons = compare_recovered_invariants_against_curve(
        curve,
        &snapshot.recovered_invariants,
        tolerance,
    )?;

    Ok(InverseUniformizationJValidationReport {
        curve: curve.clone(),
        tau: tau.clone(),
        lattice: snapshot.lattice,
        recovered_invariants: snapshot.recovered_invariants,
        comparison: comparisons.j,
    })
}

/// Validates the natural `τ` produced by a [`TauRecoveryReport`] against the
/// original curve-side `j`-invariant.
///
/// This is a convenience wrapper around
/// [`validate_recovered_tau_by_j_invariant`] that avoids manually extracting
/// the curve and the recovered natural `τ`.
pub fn validate_tau_recovery_report_by_j_invariant(
    report: &TauRecoveryReport,
    lattice_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
    let tau = report.tau();
    validate_recovered_tau_by_j_invariant(report.curve(), &tau, lattice_truncation, tolerance)
}

/// Validates the canonicalized `τ` produced by a
/// [`CanonicalTauRecoveryReport`] against the original curve-side
/// `j`-invariant.
///
/// Since `j` is modular-invariant, this checks the same modular isomorphism
/// class after the additional `SL₂(ℤ)` reduction step that sends the natural
/// recovered `τ` to the standard fundamental domain.
pub fn validate_canonical_tau_recovery_by_j_invariant(
    report: &CanonicalTauRecoveryReport,
    lattice_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
    validate_recovered_tau_by_j_invariant(
        report.curve(),
        report.canonical_tau(),
        lattice_truncation,
        tolerance,
    )
}
