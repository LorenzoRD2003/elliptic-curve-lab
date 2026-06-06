use crate::elliptic_curves::analytic::lattice::HasAnalyticLatticeContext;
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance,
    ComplexLattice, LatticeSumTruncation, RecoveredPeriodBasis, UpperHalfPlanePoint,
};
use crate::numerics::ComplexApproxComparison;

use crate::elliptic_curves::analytic::inverse_uniformization::validation_shared::{
    classify_invariant_recovery_interpretation, compare_recovered_invariants_against_curve,
    recover_invariant_snapshot_from_lattice,
};

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

    /// Returns the direct comparison between recovered and curve-side `g₂`.
    pub(crate) fn g2_comparison(&self) -> &ComplexApproxComparison {
        &self.g2_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `g₃`.
    pub(crate) fn g3_comparison(&self) -> &ComplexApproxComparison {
        &self.g3_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `Δ`.
    pub(crate) fn discriminant_comparison(&self) -> &ComplexApproxComparison {
        &self.discriminant_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `j`.
    pub(crate) fn j_comparison(&self) -> &ComplexApproxComparison {
        &self.j_comparison
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
        self.recovered_invariants.truncation
    }

    /// Returns the shared comparison tolerance.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.j_comparison.tolerance()
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

/// Validates a recovered lattice directly at the level of analytic invariants.
///
/// Given one recovered period basis `Λ_rec = ℤω₁ + ℤω₂`, this routine
/// recomputes
///
/// - `g₂(Λ_rec)`
/// - `g₃(Λ_rec)`
/// - `Δ(Λ_rec) = g₂(Λ_rec)^3 - 27 g₃(Λ_rec)^2`
/// - `j(Λ_rec)`
///
/// by finite lattice sums, and compares them against the target curve-side
/// invariants. This is a stricter validation than checking `j` alone, because
/// `g₂`, `g₃`, and `Δ` depend on the overall lattice scale.
///
/// If the recovered basis differs from the target analytic normalization by a
/// homothety `Λ_rec = αΛ`, then one should expect
///
/// - `g₂(Λ_rec) = α^{-4} g₂(Λ)`
/// - `g₃(Λ_rec) = α^{-6} g₃(Λ)`
/// - `Δ(Λ_rec) = α^{-12} Δ(Λ)`
/// - `j(Λ_rec) = j(Λ)`
///
/// so the report distinguishes:
///
/// - full direct agreement
/// - agreement only at the modular-invariant level through `j`
/// - genuine inconsistency even at the `j` level
///
/// Complexity: `Θ(r²)` in the lattice truncation radius `r`, since the
/// dominant work is the recomputation of the truncated Eisenstein sums on the
/// recovered lattice.
pub fn validate_recovered_lattice_invariants(
    curve: &AnalyticWeierstrassCurve,
    periods: &RecoveredPeriodBasis,
    lattice_truncation: LatticeSumTruncation,
    tolerance: ApproxTolerance,
) -> Result<InvariantRecoveryValidationReport, AnalyticCurveError> {
    let snapshot =
        recover_invariant_snapshot_from_lattice(periods.lattice().clone(), lattice_truncation)?;
    let comparisons = compare_recovered_invariants_against_curve(
        curve,
        &snapshot.recovered_invariants,
        tolerance,
    )?;

    Ok(InvariantRecoveryValidationReport {
        curve: curve.clone(),
        periods: periods.clone(),
        tau: periods.tau(),
        recovered_invariants: snapshot.recovered_invariants,
        interpretation: classify_invariant_recovery_interpretation(&comparisons),
        g2_comparison: comparisons.g2,
        g3_comparison: comparisons.g3,
        discriminant_comparison: comparisons.discriminant,
        j_comparison: comparisons.j,
    })
}
