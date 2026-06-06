use num_complex::Complex64;

use super::super::AnalyticCurveError;
use super::{
    ComplexAgmConfig, ComplexAgmStatus, LegendreParameter, LegendreReduction, PeriodRecoveryConfig,
    PeriodRecoveryMethod, PeriodRecoveryStatus, complex_agm,
};
use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Structured numerical metadata for one complete-elliptic-integral
/// approximation of the first kind.
#[derive(Clone, Debug, PartialEq)]
pub struct CompleteEllipticIntegralKMetadata {
    resolved_method: PeriodRecoveryMethod,
    status: PeriodRecoveryStatus,
    tolerance: ApproxTolerance,
    agm_iterations_used: usize,
    complementary_square_root: Complex64,
    used_principal_complementary_branch: bool,
}

impl CompleteEllipticIntegralKMetadata {
    /// Builds one complete-elliptic-integral metadata bundle from explicit
    /// fields.
    pub fn new(
        resolved_method: PeriodRecoveryMethod,
        status: PeriodRecoveryStatus,
        tolerance: ApproxTolerance,
        agm_iterations_used: usize,
        complementary_square_root: Complex64,
        used_principal_complementary_branch: bool,
    ) -> Self {
        Self {
            resolved_method,
            status,
            tolerance,
            agm_iterations_used,
            complementary_square_root,
            used_principal_complementary_branch,
        }
    }

    pub fn resolved_method(&self) -> PeriodRecoveryMethod {
        self.resolved_method
    }

    pub fn status(&self) -> PeriodRecoveryStatus {
        self.status
    }

    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    pub fn agm_iterations_used(&self) -> usize {
        self.agm_iterations_used
    }

    pub fn complementary_square_root(&self) -> &Complex64 {
        &self.complementary_square_root
    }

    pub fn used_principal_complementary_branch(&self) -> bool {
        self.used_principal_complementary_branch
    }

    pub fn succeeded(&self) -> bool {
        self.status == PeriodRecoveryStatus::Succeeded
    }
}

/// Approximation bundle for the complete elliptic integral of the first kind.
#[derive(Clone, Debug, PartialEq)]
pub struct CompleteEllipticIntegralKApprox {
    parameter: LegendreParameter,
    value: Complex64,
    metadata: CompleteEllipticIntegralKMetadata,
}

impl CompleteEllipticIntegralKApprox {
    pub fn new(
        parameter: LegendreParameter,
        value: Complex64,
        metadata: CompleteEllipticIntegralKMetadata,
    ) -> Self {
        Self {
            parameter,
            value,
            metadata,
        }
    }

    pub fn parameter(&self) -> &LegendreParameter {
        &self.parameter
    }

    pub fn value(&self) -> &Complex64 {
        &self.value
    }

    pub fn metadata(&self) -> &CompleteEllipticIntegralKMetadata {
        &self.metadata
    }
}

/// One high-level Legendre-side period report assembled from complete elliptic
/// integrals of the first kind.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendrePeriodIntegralReport {
    pub lambda: LegendreParameter,
    pub k_lambda: CompleteEllipticIntegralKApprox,
    pub k_complementary: CompleteEllipticIntegralKApprox,
    pub tau_candidate: Complex64,
}

/// Approximates the complete elliptic integral of the first kind from the
/// raw parameter `m = k²`.
pub fn complete_elliptic_integral_k_from_m(
    m: Complex64,
    config: ComplexAgmConfig,
) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
    let parameter = validate_complete_elliptic_integral_parameter(m)?;
    let (complementary_square_root, used_principal_complementary_branch) =
        select_complementary_square_root(m)?;

    evaluate_complete_elliptic_integral_via_agm(
        parameter,
        complementary_square_root,
        used_principal_complementary_branch,
        config,
    )
}

/// Approximates the complete elliptic integral of the first kind from a
/// validated Legendre parameter `λ`.
pub fn complete_elliptic_integral_k_from_lambda(
    parameter: &LegendreParameter,
    config: ComplexAgmConfig,
) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
    complete_elliptic_integral_k_from_m(*parameter.lambda(), config)
}

/// Approximates the complementary complete elliptic integral of the first
/// kind from the raw parameter `m = k²`.
pub fn complementary_complete_elliptic_integral_k_from_m(
    m: Complex64,
    config: ComplexAgmConfig,
) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
    complete_elliptic_integral_k_from_m(Complex64::new(1.0, 0.0) - m, config)
}

/// Approximates the complementary complete elliptic integral of the first
/// kind from a validated Legendre parameter `λ`.
pub fn complementary_complete_elliptic_integral_k_from_lambda(
    parameter: &LegendreParameter,
    config: ComplexAgmConfig,
) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
    complementary_complete_elliptic_integral_k_from_m(*parameter.lambda(), config)
}

/// Builds a Legendre-side period report from one explicit Legendre reduction.
///
/// The current period-ratio candidate is `τ_candidate = i K(1 - λ) / K(λ)`.
pub fn legendre_period_integral_report(
    reduction: &LegendreReduction,
    config: PeriodRecoveryConfig,
) -> Result<LegendrePeriodIntegralReport, AnalyticCurveError> {
    let agm_config = ComplexAgmConfig::from_period_recovery_config(config);
    let lambda = reduction.parameter().clone();
    let k_lambda = complete_elliptic_integral_k_from_lambda(&lambda, agm_config)?;
    let k_complementary =
        complementary_complete_elliptic_integral_k_from_lambda(&lambda, agm_config)?;

    if ComplexApprox::is_zero_with_tolerance(k_lambda.value(), agm_config.tolerance()) {
        return Err(AnalyticCurveError::PeriodRecoveryFailed);
    }

    let tau_candidate = Complex64::new(0.0, 1.0) * (*k_complementary.value() / *k_lambda.value());

    if !tau_candidate.is_finite() {
        return Err(AnalyticCurveError::PeriodRecoveryFailed);
    }

    Ok(LegendrePeriodIntegralReport {
        lambda,
        k_lambda,
        k_complementary,
        tau_candidate,
    })
}

fn validate_complete_elliptic_integral_parameter(
    m: Complex64,
) -> Result<LegendreParameter, AnalyticCurveError> {
    if !m.is_finite() {
        return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
    }

    LegendreParameter::new(m).map_err(|_| AnalyticCurveError::InvalidEllipticIntegralInput)
}

fn select_complementary_square_root(m: Complex64) -> Result<(Complex64, bool), AnalyticCurveError> {
    let principal = (Complex64::new(1.0, 0.0) - m).sqrt();
    if !principal.is_finite() {
        return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
    }

    let negated = -principal;
    let one = Complex64::new(1.0, 0.0);
    let principal_gap = (one - principal).norm();
    let negated_gap = (one - negated).norm();

    if principal_gap <= negated_gap {
        Ok((principal, true))
    } else {
        Ok((negated, false))
    }
}

fn evaluate_complete_elliptic_integral_via_agm(
    parameter: LegendreParameter,
    complementary_square_root: Complex64,
    used_principal_complementary_branch: bool,
    config: ComplexAgmConfig,
) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
    let agm_result = complex_agm(Complex64::new(1.0, 0.0), complementary_square_root, config)?;

    if ComplexApprox::is_zero_with_tolerance(agm_result.agm(), config.tolerance()) {
        return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
    }

    let value =
        Complex64::new(std::f64::consts::PI, 0.0) / (Complex64::new(2.0, 0.0) * *agm_result.agm());

    if !value.is_finite() {
        return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
    }

    let metadata = build_complete_elliptic_integral_k_metadata(
        &agm_result,
        config.tolerance(),
        complementary_square_root,
        used_principal_complementary_branch,
    );

    Ok(CompleteEllipticIntegralKApprox::new(
        parameter, value, metadata,
    ))
}

fn build_complete_elliptic_integral_k_metadata(
    agm_result: &super::ComplexAgmResult,
    tolerance: ApproxTolerance,
    complementary_square_root: Complex64,
    used_principal_complementary_branch: bool,
) -> CompleteEllipticIntegralKMetadata {
    let status = match agm_result.status() {
        ComplexAgmStatus::Succeeded => PeriodRecoveryStatus::Succeeded,
        ComplexAgmStatus::HitIterationLimit => PeriodRecoveryStatus::HitIterationLimit,
    };

    CompleteEllipticIntegralKMetadata::new(
        PeriodRecoveryMethod::AgmViaLegendre,
        status,
        tolerance,
        agm_result.iterations_used(),
        complementary_square_root,
        used_principal_complementary_branch,
    )
}
