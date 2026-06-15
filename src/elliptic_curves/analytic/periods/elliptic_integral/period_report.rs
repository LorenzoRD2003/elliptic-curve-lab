use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    periods::{
        CompleteEllipticIntegralKApprox, LegendreParameter, LegendreReduction,
        PeriodRecoveryConfig, elliptic_integral::ComplexAgmConfig,
    },
};
use crate::fields::complex_approx::ComplexApprox;

/// One high-level Legendre-side period report assembled from complete elliptic
/// integrals of the first kind.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendrePeriodIntegralReport {
    lambda: LegendreParameter,
    k_lambda: CompleteEllipticIntegralKApprox,
    k_complementary: CompleteEllipticIntegralKApprox,
    tau_candidate: Complex64,
}

impl LegendrePeriodIntegralReport {
    pub(crate) fn new(
        lambda: LegendreParameter,
        k_lambda: CompleteEllipticIntegralKApprox,
        k_complementary: CompleteEllipticIntegralKApprox,
        tau_candidate: Complex64,
    ) -> Self {
        Self {
            lambda,
            k_lambda,
            k_complementary,
            tau_candidate,
        }
    }

    pub fn lambda(&self) -> &LegendreParameter {
        &self.lambda
    }

    pub fn k_lambda(&self) -> &CompleteEllipticIntegralKApprox {
        &self.k_lambda
    }

    pub fn k_complementary(&self) -> &CompleteEllipticIntegralKApprox {
        &self.k_complementary
    }

    pub fn tau_candidate(&self) -> &Complex64 {
        &self.tau_candidate
    }
}

impl LegendreReduction {
    /// Builds a Legendre-side period report from this explicit reduction.
    pub fn period_integral_report(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<LegendrePeriodIntegralReport, AnalyticCurveError> {
        let agm_config = ComplexAgmConfig::from_period_recovery_config(config);
        let lambda = self.parameter().clone();
        let k_lambda = lambda.complete_elliptic_integral_k(agm_config)?;
        let k_complementary = lambda.complementary_complete_elliptic_integral_k(agm_config)?;

        if ComplexApprox::is_zero_with_tolerance(k_lambda.value(), agm_config.tolerance()) {
            return Err(AnalyticCurveError::PeriodRecoveryFailed);
        }

        let tau_candidate =
            Complex64::new(0.0, 1.0) * (*k_complementary.value() / *k_lambda.value());

        if !tau_candidate.is_finite() {
            return Err(AnalyticCurveError::PeriodRecoveryFailed);
        }

        Ok(LegendrePeriodIntegralReport::new(
            lambda,
            k_lambda,
            k_complementary,
            tau_candidate,
        ))
    }
}
