use crate::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve, FundamentalDomainReductionReport, NumericalRecoveryMetadata,
    RecoveredPeriodBasisReport, TauRecoveryReport, UpperHalfPlanePoint, WeierstrassCubicRoots,
    periods::RecoveredPeriodBasis,
};

/// A `τ`-recovery report together with an explicit canonicalization to the
/// standard fundamental domain of `SL₂(ℤ)`.
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalTauRecoveryReport {
    tau_recovery_report: TauRecoveryReport,
    fundamental_domain_reduction: FundamentalDomainReductionReport,
}

impl CanonicalTauRecoveryReport {
    pub(crate) fn new(
        tau_recovery_report: TauRecoveryReport,
        fundamental_domain_reduction: FundamentalDomainReductionReport,
    ) -> Self {
        Self {
            tau_recovery_report,
            fundamental_domain_reduction,
        }
    }

    pub fn tau_recovery_report(&self) -> &TauRecoveryReport {
        &self.tau_recovery_report
    }

    pub fn fundamental_domain_reduction(&self) -> &FundamentalDomainReductionReport {
        &self.fundamental_domain_reduction
    }

    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        self.tau_recovery_report.curve()
    }

    pub fn roots(&self) -> &WeierstrassCubicRoots {
        self.tau_recovery_report.roots()
    }

    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        self.tau_recovery_report.basis_report()
    }

    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.tau_recovery_report.periods()
    }

    pub fn original_tau(&self) -> UpperHalfPlanePoint {
        self.tau_recovery_report.tau()
    }

    pub fn canonical_tau(&self) -> &UpperHalfPlanePoint {
        self.fundamental_domain_reduction.reduced_tau()
    }

    pub fn accumulated_matrix(&self) -> crate::elliptic_curves::analytic::ModularMatrix {
        self.fundamental_domain_reduction.accumulated_matrix()
    }

    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        self.tau_recovery_report.metadata()
    }
}
