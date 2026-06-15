use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, LatticeSumTruncation,
    RecoveredPeriodBasis,
    inverse_uniformization::{
        lattice_invariants::InvariantRecoveryValidationReport,
        validation_shared::{
            classify_invariant_recovery_interpretation, compare_recovered_invariants_against_curve,
            recover_invariant_snapshot_from_lattice,
        },
    },
};

impl AnalyticWeierstrassCurve {
    /// Validates a recovered period basis against this analytic curve through
    /// the scale-sensitive invariants `g₂`, `g₃`, `Δ`, and `j`.
    pub fn validate_recovered_lattice_invariants(
        &self,
        periods: &RecoveredPeriodBasis,
        lattice_truncation: LatticeSumTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<InvariantRecoveryValidationReport, AnalyticCurveError> {
        let snapshot =
            recover_invariant_snapshot_from_lattice(periods.lattice().clone(), lattice_truncation)?;
        let comparisons = compare_recovered_invariants_against_curve(
            self,
            &snapshot.recovered_invariants,
            tolerance,
        )?;
        let interpretation = classify_invariant_recovery_interpretation(&comparisons);

        Ok(InvariantRecoveryValidationReport::new(
            self.clone(),
            periods.clone(),
            periods.tau(),
            snapshot.recovered_invariants,
            comparisons.g2,
            comparisons.g3,
            comparisons.discriminant,
            comparisons.j,
            interpretation,
        ))
    }
}
