use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ApproxTolerance, CanonicalTauRecoveryReport,
    LatticeSumTruncation, TauRecoveryReport, UpperHalfPlanePoint,
    inverse_uniformization::{
        j_validation::InverseUniformizationJValidationReport,
        validation_shared::{
            compare_recovered_invariants_against_curve, recover_invariant_snapshot_from_tau,
        },
    },
};

impl AnalyticWeierstrassCurve {
    /// Validates one recovered upper-half-plane parameter `τ` against this
    /// analytic curve by comparing modular `j`-invariants.
    pub fn validate_recovered_tau_by_j_invariant(
        &self,
        tau: &UpperHalfPlanePoint,
        lattice_truncation: LatticeSumTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
        let snapshot = recover_invariant_snapshot_from_tau(tau, lattice_truncation)?;
        let comparisons = compare_recovered_invariants_against_curve(
            self,
            &snapshot.recovered_invariants,
            tolerance,
        )?;

        Ok(InverseUniformizationJValidationReport::new(
            self.clone(),
            tau.clone(),
            snapshot.lattice,
            snapshot.recovered_invariants,
            comparisons.j,
        ))
    }
}

impl TauRecoveryReport {
    /// Validates the naturally recovered `τ` against the original curve-side
    /// `j`-invariant.
    pub fn validate_by_j_invariant(
        &self,
        lattice_truncation: LatticeSumTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
        let tau = self.tau();
        self.curve()
            .validate_recovered_tau_by_j_invariant(&tau, lattice_truncation, tolerance)
    }
}

impl CanonicalTauRecoveryReport {
    /// Validates the canonicalized recovered `τ` against the original
    /// curve-side `j`-invariant.
    pub fn validate_by_j_invariant(
        &self,
        lattice_truncation: LatticeSumTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<InverseUniformizationJValidationReport, AnalyticCurveError> {
        self.curve().validate_recovered_tau_by_j_invariant(
            self.canonical_tau(),
            lattice_truncation,
            tolerance,
        )
    }
}
