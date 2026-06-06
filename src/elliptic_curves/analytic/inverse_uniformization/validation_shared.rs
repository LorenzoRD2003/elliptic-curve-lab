use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ApproxTolerance,
    ComplexLattice, InvariantRecoveryInterpretation, LatticeSumTruncation, UpperHalfPlanePoint,
    analytic_discriminant, analytic_invariants,
};
use crate::numerics::ComplexApproxComparison;

/// Cached lattice plus recomputed invariant bundle shared by the current
/// inverse-uniformization validation helpers.
///
/// Both the `j`-only validation and the full scale-sensitive invariant
/// validation start from the same numerical subproblem:
///
/// 1. build a lattice,
/// 2. recompute `g₂`, `g₃`, `Δ`, and `j` by truncated lattice sums,
/// 3. compare those values against the target curve-side invariants.
///
/// Keeping this helper private but shared avoids duplicating the same
/// truncation-side work in two sibling modules.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RecoveredLatticeInvariantSnapshot {
    pub(crate) lattice: ComplexLattice,
    pub(crate) recovered_invariants: AnalyticInvariants,
}

/// Direct comparisons between recomputed lattice-side invariants and the
/// target curve-side invariants.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CurveInvariantComparisons {
    pub(crate) g2: ComplexApproxComparison,
    pub(crate) g3: ComplexApproxComparison,
    pub(crate) discriminant: ComplexApproxComparison,
    pub(crate) j: ComplexApproxComparison,
}

/// Recomputes the analytic invariants attached to the standard lattice
/// `Λ_τ = ℤ + ℤτ`.
pub(crate) fn recover_invariant_snapshot_from_tau(
    tau: &UpperHalfPlanePoint,
    truncation: LatticeSumTruncation,
) -> Result<RecoveredLatticeInvariantSnapshot, AnalyticCurveError> {
    let lattice = ComplexLattice::from_tau(tau.clone());
    recover_invariant_snapshot_from_lattice(lattice, truncation)
}

/// Recomputes the analytic invariants attached to an explicit lattice.
pub(crate) fn recover_invariant_snapshot_from_lattice(
    lattice: ComplexLattice,
    truncation: LatticeSumTruncation,
) -> Result<RecoveredLatticeInvariantSnapshot, AnalyticCurveError> {
    let recovered_invariants = analytic_invariants(&lattice, truncation)?;

    Ok(RecoveredLatticeInvariantSnapshot {
        lattice,
        recovered_invariants,
    })
}

/// Compares recomputed lattice-side invariants against one target analytic
/// Weierstrass curve.
pub(crate) fn compare_recovered_invariants_against_curve(
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

/// Classifies whether the recovered lattice agreed directly, only up to
/// modular class, or not even at the `j`-invariant level.
pub(crate) fn classify_invariant_recovery_interpretation(
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
