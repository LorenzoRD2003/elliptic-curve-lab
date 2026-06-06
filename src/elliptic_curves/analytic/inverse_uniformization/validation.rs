use num_complex::Complex64;

use super::super::periods::{CanonicalTauRecoveryReport, RecoveredPeriodBasis, TauRecoveryReport};
use super::super::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance,
    ComplexApproxComparison, ComplexLattice, HasAnalyticLatticeContext, HasComplexApproxComparison,
    LatticeSumTruncation, UpperHalfPlanePoint, analytic_discriminant, analytic_invariants,
};

#[derive(Clone, Debug, PartialEq)]
struct RecoveredLatticeInvariantSnapshot {
    lattice: ComplexLattice,
    recovered_invariants: AnalyticInvariants,
}

#[derive(Clone, Debug, PartialEq)]
struct CurveInvariantComparisons {
    g2: ComplexApproxComparison,
    g3: ComplexApproxComparison,
    discriminant: ComplexApproxComparison,
    j: ComplexApproxComparison,
}

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
    pub fn g2_comparison(&self) -> &ComplexApproxComparison {
        &self.g2_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `g₃`.
    pub fn g3_comparison(&self) -> &ComplexApproxComparison {
        &self.g3_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `Δ`.
    pub fn discriminant_comparison(&self) -> &ComplexApproxComparison {
        &self.discriminant_comparison
    }

    /// Returns the direct comparison between recovered and curve-side `j`.
    pub fn j_comparison(&self) -> &ComplexApproxComparison {
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

fn recover_invariant_snapshot_from_tau(
    tau: &UpperHalfPlanePoint,
    truncation: LatticeSumTruncation,
) -> Result<RecoveredLatticeInvariantSnapshot, AnalyticCurveError> {
    let lattice = ComplexLattice::from_tau(tau.clone());
    recover_invariant_snapshot_from_lattice(lattice, truncation)
}

fn recover_invariant_snapshot_from_lattice(
    lattice: ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<RecoveredLatticeInvariantSnapshot, AnalyticCurveError> {
    let recovered_invariants = analytic_invariants(&lattice, truncation)?;

    Ok(RecoveredLatticeInvariantSnapshot {
        lattice,
        recovered_invariants,
    })
}

fn compare_recovered_invariants_against_curve(
    curve: &AnalyticWeierstrassCurve,
    recovered_invariants: &AnalyticInvariants,
    tolerance: ApproxTolerance,
) -> Result<CurveInvariantComparisons, AnalyticCurveError> {
    let curve_discriminant = analytic_discriminant(curve.g2(), curve.g3());
    let curve_j = curve.j_invariant()?;

    Ok(CurveInvariantComparisons {
        g2: ComplexApproxComparison::new(recovered_invariants.g2, *curve.g2(), tolerance),
        g3: ComplexApproxComparison::new(recovered_invariants.g3, *curve.g3(), tolerance),
        discriminant: ComplexApproxComparison::new(
            recovered_invariants.discriminant,
            curve_discriminant,
            tolerance,
        ),
        j: ComplexApproxComparison::new(recovered_invariants.j_invariant, curve_j, tolerance),
    })
}

fn classify_invariant_recovery_interpretation(
    comparisons: &CurveInvariantComparisons,
) -> InvariantRecoveryInterpretation {
    let direct_scale_sensitive_agreement = comparisons.g2.agrees_approximately()
        && comparisons.g3.agrees_approximately()
        && comparisons.discriminant.agrees_approximately();

    if direct_scale_sensitive_agreement && comparisons.j.agrees_approximately() {
        InvariantRecoveryInterpretation::DirectAgreement
    } else if comparisons.j.agrees_approximately() {
        InvariantRecoveryInterpretation::SameModularClassButScaleSensitiveMismatch
    } else {
        InvariantRecoveryInterpretation::Inconsistent
    }
}
