use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, LatticeSumTruncation, UpperHalfPlanePoint,
    analytic_invariants_from_tau,
};
use crate::fields::ComplexApprox;

use super::{JInvariantQExpansion, QExpansionTruncation};

/// Structured comparison between the Eisenstein-sum and `q`-expansion routes
/// to the analytic modular `j`-invariant.
#[derive(Clone, Debug, PartialEq)]
pub struct JInvariantComparisonReport {
    tau: UpperHalfPlanePoint,
    eisenstein_j: Complex64,
    q_expansion_j: Complex64,
    difference: Complex64,
    agrees_approximately: bool,
    lattice_truncation: LatticeSumTruncation,
    q_truncation: QExpansionTruncation,
    tolerance: ApproxTolerance,
}

impl JInvariantComparisonReport {
    /// Returns the upper-half-plane parameter `τ` used in both computations.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the `j`-invariant obtained from truncated Eisenstein sums.
    pub fn eisenstein_j(&self) -> &Complex64 {
        &self.eisenstein_j
    }

    /// Returns the `j`-invariant obtained from the truncated `q`-expansion.
    pub fn q_expansion_j(&self) -> &Complex64 {
        &self.q_expansion_j
    }

    /// Returns the residual `j_Eisenstein - j_q`.
    pub fn difference(&self) -> &Complex64 {
        &self.difference
    }

    /// Returns the Euclidean norm `|j_Eisenstein - j_q|`.
    pub fn absolute_difference(&self) -> f64 {
        self.difference.norm()
    }

    /// Returns whether the two approximations agreed under the supplied
    /// tolerance policy.
    pub fn agrees_approximately(&self) -> bool {
        self.agrees_approximately
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
        self.tolerance
    }
}

/// Compares the two current analytic routes to the modular `j`-invariant:
///
/// - truncated Eisenstein sums on the lattice `Λ_τ = ℤ + ℤτ`
/// - truncated cusp expansion in `q = e^{2π i τ}`
///
/// This is an educational numerical experiment rather than a certified
/// modular-forms routine. Its quality depends both on the lattice truncation
/// radius and on how small `|q|` is for the chosen `τ`.
///
/// Complexity: `Θ(r² + N)`, where `r` is the lattice truncation radius and
/// `N = q_truncation.terms()`.
pub fn compare_j_from_eisenstein_and_q_expansion(
    tau: UpperHalfPlanePoint,
    lattice_truncation: LatticeSumTruncation,
    q_truncation: QExpansionTruncation,
    tolerance: ApproxTolerance,
) -> Result<JInvariantComparisonReport, AnalyticCurveError> {
    let invariants = analytic_invariants_from_tau(&tau, lattice_truncation)?;
    let q_approximation = JInvariantQExpansion::from_tau(tau.clone(), q_truncation)?;
    let comparison = ComplexApprox::comparison_report(
        &invariants.j_invariant,
        q_approximation.value(),
        tolerance,
    );

    Ok(JInvariantComparisonReport {
        tau,
        eisenstein_j: invariants.j_invariant,
        q_expansion_j: *q_approximation.value(),
        difference: comparison.difference,
        agrees_approximately: comparison.is_close,
        lattice_truncation,
        q_truncation,
        tolerance,
    })
}
