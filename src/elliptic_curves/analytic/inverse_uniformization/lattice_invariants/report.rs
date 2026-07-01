use crate::elliptic_curves::analytic::{
    AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
    LatticeSumTruncation, RecoveredPeriodBasis, UpperHalfPlanePoint,
    lattice::HasAnalyticLatticeContext,
};
use crate::numerics::ComplexApproxComparison;

/// Interprets how a recovered lattice compares to the target analytic curve.
///
/// The key mathematical distinction is that `j` is modular- and
/// homothety-invariant, while `g₂`, `g₃`, and `Δ` are scale-sensitive:
///
/// - if `Λ' = αΛ`, then
///   `g₂(Λ') = α^{-4} g₂(Λ)`,
///   `g₃(Λ') = α^{-6} g₃(Λ)`,
///   `Δ(Λ') = α^{-12} Δ(Λ)`,
///   but `j(Λ') = j(Λ)`.
///
/// So agreement of `j` alone can still mean that the recovered lattice is in
/// the correct modular class while using a different overall homothety.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InvariantRecoveryInterpretation {
    /// The recovered lattice agrees directly with the target curve at the
    /// level of `g₂`, `g₃`, `Δ`, and `j`.
    DirectAgreement,
    /// The recovered lattice matches the target curve only at the modular
    /// level through `j`; at least one of `g₂`, `g₃`, or `Δ` still differs.
    SameModularClassButScaleSensitiveMismatch,
    /// Even the modular-invariant quantity `j` failed to agree numerically.
    Inconsistent,
}

/// Validation report that compares the scale-sensitive analytic invariants
/// recovered from one lattice against the target curve-side invariants.
///
/// Unlike a pure `j`-comparison, this report is deliberately sensitive to the
/// chosen homothety normalization of the recovered lattice. It therefore keeps
/// both stories visible:
///
/// - direct comparison of `g₂`, `g₃`, `Δ`, and `j`
/// - an interpretation explaining whether any mismatch can still be
///   understood as “same modular class, different overall scale”
#[derive(Clone, Debug, PartialEq)]
pub struct InvariantRecoveryValidationReport {
    curve: AnalyticWeierstrassCurve,
    periods: RecoveredPeriodBasis,
    tau: UpperHalfPlanePoint,
    recovered_invariants: AnalyticInvariants,
    g2_comparison: ComplexApproxComparison,
    g3_comparison: ComplexApproxComparison,
    discriminant_comparison: ComplexApproxComparison,
    j_comparison: ComplexApproxComparison,
    interpretation: InvariantRecoveryInterpretation,
}

impl InvariantRecoveryValidationReport {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        curve: AnalyticWeierstrassCurve,
        periods: RecoveredPeriodBasis,
        tau: UpperHalfPlanePoint,
        recovered_invariants: AnalyticInvariants,
        g2_comparison: ComplexApproxComparison,
        g3_comparison: ComplexApproxComparison,
        discriminant_comparison: ComplexApproxComparison,
        j_comparison: ComplexApproxComparison,
        interpretation: InvariantRecoveryInterpretation,
    ) -> Self {
        Self {
            curve,
            periods,
            tau,
            recovered_invariants,
            g2_comparison,
            g3_comparison,
            discriminant_comparison,
            j_comparison,
            interpretation,
        }
    }

    /// Returns the target analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the recovered period basis whose lattice is being validated.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        &self.periods
    }

    /// Returns the corresponding period ratio `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the recovered lattice invariants.
    pub fn recovered_invariants(&self) -> &AnalyticInvariants {
        &self.recovered_invariants
    }

    /// Returns the interpretation of the combined comparison outcome.
    pub fn interpretation(&self) -> InvariantRecoveryInterpretation {
        self.interpretation
    }

    /// Returns whether the scale-sensitive invariants `g₂`, `g₃`, and `Δ`
    /// all agreed directly.
    pub fn direct_scale_sensitive_agreement(&self) -> bool {
        self.g2_comparison.agrees_approximately()
            && self.g3_comparison.agrees_approximately()
            && self.discriminant_comparison.agrees_approximately()
    }

    /// Returns whether the modular-invariant quantity `j` agreed.
    pub fn same_j_invariant_approximately(&self) -> bool {
        self.j_comparison.agrees_approximately()
    }

    /// Returns the lattice-sum truncation used on the recovered lattice side.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.recovered_invariants.truncation()
    }

    /// Returns the shared comparison tolerance.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.j_comparison.tolerance()
    }

    /// Returns the direct comparison between recovered and curve-side `g₂`.
    #[cfg(feature = "visualization")]
    pub(crate) fn g2_comparison(&self) -> &ComplexApproxComparison {
        &self.g2_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `g₃`.
    #[cfg(feature = "visualization")]
    pub(crate) fn g3_comparison(&self) -> &ComplexApproxComparison {
        &self.g3_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `Δ`.
    #[cfg(feature = "visualization")]
    pub(crate) fn discriminant_comparison(&self) -> &ComplexApproxComparison {
        &self.discriminant_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `j`.
    #[cfg(feature = "visualization")]
    pub(crate) fn j_comparison(&self) -> &ComplexApproxComparison {
        &self.j_comparison
    }
}

impl HasAnalyticLatticeContext for InvariantRecoveryValidationReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        self.periods.lattice()
    }
}
