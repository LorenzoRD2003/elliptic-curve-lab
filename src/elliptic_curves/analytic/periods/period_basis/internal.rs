use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, LegendrePeriodIntegralReport, LegendreReduction, PeriodRecoveryConfig,
    periods::RecoveredPeriodBasis,
};
use crate::fields::complex_approx::ComplexApprox;

pub(super) fn build_period_basis_from_legendre_data(
    reduction: &LegendreReduction,
    integral_report: &LegendrePeriodIntegralReport,
    config: PeriodRecoveryConfig,
) -> Result<RecoveredPeriodBasis, AnalyticCurveError> {
    let differential_scale = reduction.invariant_differential_scale();
    let omega1 =
        Complex64::new(4.0, 0.0) * differential_scale * *integral_report.k_lambda().value();
    let omega2 =
        Complex64::new(0.0, 4.0) * differential_scale * *integral_report.k_complementary().value();

    if !omega1.is_finite() || !omega2.is_finite() {
        return Err(AnalyticCurveError::PeriodRecoveryFailed);
    }

    let basis = RecoveredPeriodBasis::new(omega1, omega2).map_err(map_basis_construction_error)?;
    validate_period_basis_against_tau_candidate(&basis, integral_report, config)?;

    Ok(basis)
}

fn validate_period_basis_against_tau_candidate(
    basis: &RecoveredPeriodBasis,
    integral_report: &LegendrePeriodIntegralReport,
    config: PeriodRecoveryConfig,
) -> Result<(), AnalyticCurveError> {
    let tau = basis.tau();
    let tau_value = *tau.tau();
    if !ComplexApprox::eq_with_tolerance(
        &tau_value,
        integral_report.tau_candidate(),
        config.tolerance(),
    ) {
        return Err(AnalyticCurveError::PeriodValidationFailed);
    }

    Ok(())
}

fn map_basis_construction_error(error: AnalyticCurveError) -> AnalyticCurveError {
    match error {
        AnalyticCurveError::NonPositiveLatticeOrientation => {
            AnalyticCurveError::PeriodRatioNotInUpperHalfPlane
        }
        AnalyticCurveError::DegenerateLattice => AnalyticCurveError::PeriodRecoveryFailed,
        other => other,
    }
}
