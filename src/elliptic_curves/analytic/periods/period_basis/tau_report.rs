use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, FundamentalDomainReductionReport,
    FundamentalDomainReductionStatus, LegendrePeriodIntegralReport, LegendreReduction,
    NumericalRecoveryMetadata, PeriodBasisRecoveryReport, PeriodRecoveryConfig,
    RecoveredPeriodBasisReport, UpperHalfPlanePoint, WeierstrassCubicRoots,
    periods::RecoveredPeriodBasis,
};

/// A focused curve-level report for users who primarily want the recovered
/// modular parameter `τ`.
#[derive(Clone, Debug, PartialEq)]
pub struct TauRecoveryReport {
    period_basis_report: PeriodBasisRecoveryReport,
}

impl TauRecoveryReport {
    pub(crate) fn new(period_basis_report: PeriodBasisRecoveryReport) -> Self {
        Self {
            period_basis_report,
        }
    }

    pub fn period_basis_report(&self) -> &PeriodBasisRecoveryReport {
        &self.period_basis_report
    }

    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        self.period_basis_report.curve()
    }

    pub fn roots(&self) -> &WeierstrassCubicRoots {
        self.period_basis_report.roots()
    }

    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        self.period_basis_report.basis_report()
    }

    pub fn legendre_reduction(&self) -> &LegendreReduction {
        self.period_basis_report.legendre_reduction()
    }

    pub fn k_report(&self) -> &LegendrePeriodIntegralReport {
        self.period_basis_report.k_report()
    }

    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.period_basis_report.periods()
    }

    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.period_basis_report.tau()
    }

    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        self.period_basis_report.metadata()
    }

    pub(crate) fn reduce_to_standard_domain(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<FundamentalDomainReductionReport, AnalyticCurveError> {
        let reduction = self.tau().reduce_to_standard_fundamental_domain(
            config.fundamental_domain_reduction_max_steps(),
        )?;

        match reduction.status() {
            FundamentalDomainReductionStatus::AlreadyReduced
            | FundamentalDomainReductionStatus::Reduced => Ok(reduction),
            FundamentalDomainReductionStatus::StepLimitReached => {
                Err(AnalyticCurveError::PeriodValidationFailed)
            }
        }
    }
}
