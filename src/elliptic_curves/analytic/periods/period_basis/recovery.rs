use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, LegendreReduction, NumericalRecoveryMetadata,
    PeriodRecoveryConfig,
    periods::{
        CanonicalTauRecoveryReport, PeriodBasisRecoveryReport, RecoveredPeriodBasisReport,
        TauRecoveryReport, period_basis::internal::build_period_basis_from_legendre_data,
    },
};

impl LegendreReduction {
    /// Recovers one period-basis report from this explicit Legendre reduction.
    ///
    /// Complexity: `Θ(a)`, where `a = config.agm_max_iterations()`.
    pub fn recover_period_basis(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<RecoveredPeriodBasisReport, AnalyticCurveError> {
        let integral_report = self.period_integral_report(config)?;
        let basis = build_period_basis_from_legendre_data(self, &integral_report, config)?;

        Ok(RecoveredPeriodBasisReport::new(
            self.clone(),
            integral_report,
            basis,
        ))
    }
}

impl AnalyticWeierstrassCurve {
    /// Recovers one complete period-basis report directly from this analytic
    /// Weierstrass curve.
    ///
    /// Complexity: `Θ(n + a)`, where
    /// - `n = config.newton_max_iterations()`
    /// - `a = config.agm_max_iterations()`
    pub fn recover_period_basis(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<PeriodBasisRecoveryReport, AnalyticCurveError> {
        let root_report = self.recover_weierstrass_cubic_roots_with_report(config)?;
        let roots = root_report.roots().clone();
        let legendre_reduction = roots.legendre_reduction(config.tolerance())?;
        let basis_report = legendre_reduction.recover_period_basis(config)?;
        let metadata = NumericalRecoveryMetadata::from_root_and_integral_reports(
            root_report.metadata(),
            basis_report.integral_report(),
            *basis_report.tau().tau(),
        );

        Ok(PeriodBasisRecoveryReport::new(
            self.clone(),
            roots,
            basis_report,
            metadata,
        ))
    }

    /// Recovers the modular parameter `τ` directly from this analytic
    /// Weierstrass curve.
    pub fn recover_tau(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<TauRecoveryReport, AnalyticCurveError> {
        self.recover_period_basis(config)
            .map(TauRecoveryReport::new)
    }

    /// Recovers a canonically normalized modular parameter from this analytic
    /// Weierstrass curve.
    pub fn recover_canonical_tau(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<CanonicalTauRecoveryReport, AnalyticCurveError> {
        let tau_recovery_report = self.recover_tau(config)?;
        let fundamental_domain_reduction = tau_recovery_report.reduce_to_standard_domain(config)?;

        Ok(CanonicalTauRecoveryReport::new(
            tau_recovery_report,
            fundamental_domain_reduction,
        ))
    }
}
