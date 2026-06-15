use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve,
    periods::{
        CubicRootRecoveryReport, NumericalRecoveryMetadata, PeriodRecoveryConfig,
        PeriodRecoveryMethod, PeriodRecoveryStatus, WeierstrassCubicRoots,
        root_recovery::CardanoRootRecoveryDiagnostics,
    },
};
use crate::numerics::ComplexApproxComparison;

impl AnalyticWeierstrassCurve {
    /// Recovers the three roots of the Weierstrass cubic
    /// `4x^3 - g₂x - g₃ = 4(x - e₁)(x - e₂)(x - e₃)`.
    ///
    /// Complexity: `Θ(n)`, where `n = config.newton_max_iterations()`.
    pub fn recover_weierstrass_cubic_roots(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<WeierstrassCubicRoots, AnalyticCurveError> {
        self.recover_weierstrass_cubic_roots_internal(config)
            .map(|result| result.roots)
    }

    /// Recovers the Weierstrass cubic roots together with a structured
    /// reconstruction report.
    ///
    /// Complexity: `Θ(n)`, where `n = config.newton_max_iterations()`.
    pub fn recover_weierstrass_cubic_roots_with_report(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<CubicRootRecoveryReport, AnalyticCurveError> {
        let recovery = self.recover_weierstrass_cubic_roots_internal(config)?;
        let g2_comparison =
            ComplexApproxComparison::new(recovery.roots.g2(), *self.g2(), config.tolerance());
        let g3_comparison =
            ComplexApproxComparison::new(recovery.roots.g3(), *self.g3(), config.tolerance());
        let metadata = NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::Succeeded,
            recovery.newton_iterations_used,
            0,
            0,
            0,
            config.tolerance(),
            Some(recovery.validation_residual_norm),
        );
        let cardano_diagnostics = recovery
            .cardano_root_recovery_diagnostics
            .map(|diagnostics| {
                CardanoRootRecoveryDiagnostics::new(
                    diagnostics.cardano_discriminant,
                    diagnostics.cardano_product_residual_norm,
                    diagnostics.selected_u_branch_index,
                    diagnostics.selected_v_branch_index,
                )
            });

        Ok(CubicRootRecoveryReport::new(
            self.clone(),
            recovery.roots,
            g2_comparison,
            g3_comparison,
            metadata,
            cardano_diagnostics,
        ))
    }
}
