use num_complex::Complex64;
use std::f64::consts::PI;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    periods::{
        LegendreParameter, PeriodRecoveryMethod, PeriodRecoveryStatus,
        elliptic_integral::{ComplexAgmConfig, ComplexAgmResult, ComplexAgmStatus, complex_agm},
    },
};
use crate::fields::complex_approx::ComplexApprox;
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
    pub(crate) fn new(
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

    fn from_agm_result(
        agm_result: &ComplexAgmResult,
        tolerance: ApproxTolerance,
        complementary_square_root: Complex64,
        used_principal_complementary_branch: bool,
    ) -> Self {
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
}

/// Approximation bundle for the complete elliptic integral of the first kind.
#[derive(Clone, Debug, PartialEq)]
pub struct CompleteEllipticIntegralKApprox {
    parameter: LegendreParameter,
    value: Complex64,
    metadata: CompleteEllipticIntegralKMetadata,
}

impl CompleteEllipticIntegralKApprox {
    pub(crate) fn new(
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

    fn from_agm_route(
        parameter: LegendreParameter,
        complementary_square_root: Complex64,
        used_principal_complementary_branch: bool,
        config: ComplexAgmConfig,
    ) -> Result<Self, AnalyticCurveError> {
        let agm_result = complex_agm(Complex64::new(1.0, 0.0), complementary_square_root, config)?;

        if ComplexApprox::is_zero_with_tolerance(agm_result.agm(), config.tolerance()) {
            return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
        }

        let value = Complex64::new(PI, 0.0) / (Complex64::new(2.0, 0.0) * *agm_result.agm());

        if !value.is_finite() {
            return Err(AnalyticCurveError::InvalidEllipticIntegralInput);
        }

        let metadata = CompleteEllipticIntegralKMetadata::from_agm_result(
            &agm_result,
            config.tolerance(),
            complementary_square_root,
            used_principal_complementary_branch,
        );

        Ok(CompleteEllipticIntegralKApprox::new(
            parameter, value, metadata,
        ))
    }
}

impl LegendreParameter {
    /// Approximates the complete elliptic integral of the first kind `K(λ)`.
    pub fn complete_elliptic_integral_k(
        &self,
        config: ComplexAgmConfig,
    ) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
        let m = *self.lambda();
        let (complementary_square_root, used_principal_complementary_branch) =
            select_complementary_square_root(m)?;

        CompleteEllipticIntegralKApprox::from_agm_route(
            self.clone(),
            complementary_square_root,
            used_principal_complementary_branch,
            config,
        )
    }

    /// Approximates the complementary complete elliptic integral
    /// `K(1 - λ)`.
    pub fn complementary_complete_elliptic_integral_k(
        &self,
        config: ComplexAgmConfig,
    ) -> Result<CompleteEllipticIntegralKApprox, AnalyticCurveError> {
        let complementary_parameter = validate_complete_elliptic_integral_parameter(
            Complex64::new(1.0, 0.0) - *self.lambda(),
        )?;
        let (complementary_square_root, used_principal_complementary_branch) =
            select_complementary_square_root(*complementary_parameter.lambda())?;

        CompleteEllipticIntegralKApprox::from_agm_route(
            complementary_parameter,
            complementary_square_root,
            used_principal_complementary_branch,
            config,
        )
    }
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
