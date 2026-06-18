//! Abel-Jacobi recovery of torus classes from points on the analytic cubic.
//!
//! The guiding integral is
//! `z = ∫_x^∞ dt / √(4t³ - g₂ t - g₃)`,
//! interpreted modulo the recovered period lattice. In the current crate this
//! is organized as a Legendre-side numerical workflow:
//!
//! - recover or accept period data for the ambient curve,
//! - transport the point to `Y² = X(X - 1)(X - λ)`,
//! - choose one deterministic contour in the `X`-plane,
//! - follow the square-root branch continuously,
//! - integrate numerically and reduce the result modulo the lattice.
//!
//! The subfiles reflect those phases:
//!
//! - `config.rs` stores validated contour and validation policies.
//! - `contour.rs` chooses the Legendre-side contour geometry.
//! - `integration.rs` performs the branch-tracked numerical integration.
//! - `report.rs` packages the integral and point-recovery reports.
//! - `roundtrip_validation.rs` checks the recovered torus class by mapping
//!   forward again through `(℘, ℘′)`.
//!
//! This is still an educational approximation layer, so the public API keeps
//! the contour choice, branch choice, and validation residuals visible.
use num_complex::Complex64;

use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::analytic::periods::{
    LegendreReduction, PeriodRecoveryConfig, RecoveredPeriodBasis,
};
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, AnalyticWeierstrassCurve,
};
use crate::numerics::SimpsonQuadratureDomain;

mod api;
mod config;
mod contour;
mod integration;
mod report;
mod roundtrip_validation;

use contour::choose_legendre_contour;
use integration::{
    estimate_legendre_tail_correction, initialize_legendre_branch_state,
    integrate_ray_with_branch_tracking, integrate_segment_with_branch_tracking,
};
use roundtrip_validation::point_roundtrip_validation_report_for_representative;

pub use config::{
    AbelJacobiBudgets, AbelJacobiConfig, AbelJacobiRecoveryMetadata, AbelJacobiRecoveryStatus,
    AbelJacobiValidationPolicy, LegendreContourStrategy,
};
pub use report::{
    AbelJacobiContourReport, AbelJacobiInitialBranchChoice, AbelJacobiIntegralApprox,
    AbelJacobiIntegralDecomposition, AbelJacobiIntegralNumerics, AbelJacobiPointRecoveryReport,
    InverseUniformizationPointRecoveryReport,
};
pub use roundtrip_validation::{
    AbelJacobiRoundtripValidationReport, PointRoundTripValidationConfig,
    PointRoundTripValidationReport,
};

/// Approximates the Abel-Jacobi integral attached to one point of the
/// analytic Weierstrass curve.
///
/// The target quantity is the inverse-uniformization integral
/// `z = ∫_x^∞ dt / √(4 t^3 - g₂ t - g₃)`,
/// before any reduction modulo the period lattice.
///
/// For finite points, the current implementation:
///
/// 1. recovers a deterministic Legendre reduction
/// 2. transports the point to `Y² = X(X - 1)(X - λ)`
/// 3. integrates along one deterministic `segment + ray` contour in the
///    `X`-plane
/// 4. follows the square-root branch by continuity, starting from the sign
///    opposite to the input `y`-coordinate so that the convention
///    `z = ∫_x^∞ dt / sqrt(4 t^3 - g₂ t - g₃)` matches the local
///    uniformization parameter
/// 5. multiplies by the invariant-differential scale
///
/// Complexity: `Θ(n + s + r)`, where:
/// - `n = config.integration_steps`
/// - `s = config.segment_samples`
/// - `r = config.ray_samples`
///
/// The `s + r` term comes from contour scoring across a fixed finite list of
/// candidate angles.
pub(crate) fn approximate_abel_jacobi_integral_impl(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    config: AbelJacobiConfig,
) -> Result<AbelJacobiIntegralApprox, AnalyticCurveError> {
    let AffinePoint::Finite { x, y } = point else {
        return AbelJacobiIntegralApprox::new(
            curve.clone(),
            point.clone(),
            AbelJacobiContourReport::new(
                config.legendre_contour_strategy(),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                0.0,
                0.0,
                0.0,
                0.0,
            )?,
            Complex64::new(0.0, 0.0),
            AbelJacobiIntegralDecomposition::new(
                AbelJacobiInitialBranchChoice::Principal,
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
            )?,
            AbelJacobiIntegralNumerics::new(0, 0, config.tolerance()),
        );
    };

    if y.norm()
        <= config
            .tolerance()
            .absolute
            .max(config.tolerance().relative * y.norm().max(1.0))
    {
        return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
    }

    let roots = curve.recover_weierstrass_cubic_roots(config.derived_root_recovery_config())?;
    let reduction = LegendreReduction::from_roots(&roots, config.tolerance())?;
    let transformed_x = reduction.legendre_x_from_original_x(*x);
    let transformed_y = *y / reduction.legendre_y_scale();

    if !transformed_x.is_finite() || !transformed_y.is_finite() {
        return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
    }

    let contour = choose_legendre_contour(transformed_x, *reduction.parameter().lambda(), config)?;
    let contour_report = AbelJacobiContourReport::new(
        config.legendre_contour_strategy(),
        transformed_x,
        *contour.segment.end(),
        contour.ray.angle_radians(),
        contour.anchor_radius,
        contour.tail_length,
        contour.min_distance_to_branch_points,
    )?;
    let initial_branch = initialize_legendre_branch_state(
        transformed_x,
        transformed_y,
        *reduction.parameter().lambda(),
        config,
    )?;
    let initial_branch_choice = initial_branch.initial_branch_choice();

    let segment_steps = SimpsonQuadratureDomain::new(0.0, 1.0, config.integration_steps() / 2)
        .map_err(|_| AnalyticCurveError::AbelJacobiIntegrationFailed)?
        .normalized_subintervals();
    let ray_steps = SimpsonQuadratureDomain::new(
        0.0,
        contour.ray_compact_parameter_max,
        config.integration_steps() - segment_steps,
    )
    .map_err(|_| AnalyticCurveError::AbelJacobiIntegrationFailed)?
    .normalized_subintervals();
    let (segment_integral, segment_end_branch, segment_branch_adjustments) =
        integrate_segment_with_branch_tracking(
            &contour.segment,
            *reduction.parameter().lambda(),
            initial_branch,
            segment_steps,
            config,
        )?;
    let (ray_integral, ray_end_branch, ray_branch_adjustments) =
        integrate_ray_with_branch_tracking(
            &contour.ray,
            contour.ray_compact_parameter_max,
            *reduction.parameter().lambda(),
            segment_end_branch,
            ray_steps,
            config,
        )?;
    let tail_endpoint = contour
        .ray
        .point_at_compact_parameter(contour.ray_compact_parameter_max);
    let tail_correction = estimate_legendre_tail_correction(tail_endpoint, ray_end_branch)?;
    let legendre_integral = segment_integral + ray_integral + tail_correction;
    let integral_value = reduction.invariant_differential_scale() * legendre_integral;

    AbelJacobiIntegralApprox::new(
        curve.clone(),
        point.clone(),
        contour_report,
        integral_value,
        AbelJacobiIntegralDecomposition::new(
            initial_branch_choice,
            segment_integral,
            ray_integral,
            tail_correction,
        )?,
        AbelJacobiIntegralNumerics::new(
            segment_steps + ray_steps,
            segment_branch_adjustments + ray_branch_adjustments,
            config.tolerance(),
        ),
    )
}

/// Recovers the torus class of one point on an analytic Weierstrass curve,
/// assuming a period basis is already known.
///
/// The intended numerical algorithm is the Abel-Jacobi integral
/// `z = ∫_x^∞ dt / √(4 t^3 - g₂ t - g₃)`,
/// interpreted with a branch-continuation policy for the square root and then
/// reduced modulo the recovered lattice.
///
/// Special case convention:
/// the point at infinity should recover the identity class of `C / Λ`.
///
/// Complexity: `Θ(n + s + r)` for the Abel-Jacobi stage, plus the existing
/// lattice and elliptic-function validation costs, where:
/// - `n = config.integration_steps`
/// - `s = config.segment_samples`
/// - `r = config.ray_samples`
pub(crate) fn recover_torus_point_from_curve_point_with_periods_impl(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    periods: &RecoveredPeriodBasis,
    config: AbelJacobiConfig,
) -> Result<AbelJacobiPointRecoveryReport, AnalyticCurveError> {
    let report = recover_torus_point_from_curve_point_with_periods_and_validation_config(
        curve,
        point,
        periods,
        config,
        config.to_point_roundtrip_validation_config()?,
    )?;

    if !report.validation_report().agrees_approximately() {
        return Err(AnalyticCurveError::PeriodValidationFailed);
    }

    Ok(report)
}

fn recover_torus_point_from_curve_point_with_periods_and_validation_config(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    periods: &RecoveredPeriodBasis,
    config: AbelJacobiConfig,
    validation_config: PointRoundTripValidationConfig,
) -> Result<AbelJacobiPointRecoveryReport, AnalyticCurveError> {
    let integral_approx = approximate_abel_jacobi_integral_impl(curve, point, config)?;
    let torus_point = periods
        .lattice()
        .reduce_complex_point_to_torus_point(*integral_approx.value())?;
    let reduced_representative = periods
        .lattice()
        .point_from_fundamental_coordinates(torus_point.coordinate().clone());
    let validation_report = point_roundtrip_validation_report_for_representative(
        point,
        periods,
        reduced_representative,
        validation_config,
    )?;
    let status = if validation_report.agrees_approximately() {
        AbelJacobiRecoveryStatus::Succeeded
    } else {
        AbelJacobiRecoveryStatus::ValidationFailed
    };
    let metadata = AbelJacobiRecoveryMetadata::new(
        status,
        integral_approx.integration_steps_used(),
        integral_approx.branch_adjustments_used(),
        0,
        config.tolerance(),
        Some(validation_report.x_residual_norm()),
        Some(validation_report.y_residual_norm()),
    );

    AbelJacobiPointRecoveryReport::new(
        periods.clone(),
        integral_approx,
        torus_point,
        reduced_representative,
        validation_report,
        metadata,
    )
}

/// Recovers one torus class from the source curve point and immediately
/// validates the full roundtrip
///
/// `P -> z_P mod Λ -> (wp(z_P), wp'(z_P))`
///
/// using an already supplied period basis and an explicit forward-validation
/// policy.
pub(crate) fn validate_point_inverse_uniformization_roundtrip_with_periods_impl(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    periods: &RecoveredPeriodBasis,
    abel_config: AbelJacobiConfig,
    validation_config: PointRoundTripValidationConfig,
) -> Result<PointRoundTripValidationReport, AnalyticCurveError> {
    PointRoundTripValidationReport::new(
        recover_torus_point_from_curve_point_with_periods_and_validation_config(
            curve,
            point,
            periods,
            abel_config,
            validation_config,
        )?,
        validation_config,
    )
}

/// End-to-end convenience wrapper that first recovers periods from the curve
/// and then applies the Abel-Jacobi inverse map to the requested point.
pub(crate) fn recover_torus_point_from_curve_point_impl(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    period_config: PeriodRecoveryConfig,
    abel_jacobi_config: AbelJacobiConfig,
) -> Result<InverseUniformizationPointRecoveryReport, AnalyticCurveError> {
    let period_basis_report = curve.recover_period_basis(period_config)?;
    let point_recovery_report = recover_torus_point_from_curve_point_with_periods_impl(
        curve,
        point,
        period_basis_report.periods(),
        abel_jacobi_config,
    )?;

    InverseUniformizationPointRecoveryReport::new(period_basis_report, point_recovery_report)
}

/// End-to-end wrapper for the point-level inverse-uniformization roundtrip
/// experiment.
///
/// This first recovers one period basis from the curve and then runs
///
/// `P -> z_P mod Λ -> (wp(z_P), wp'(z_P))`
///
/// under the supplied forward-validation policy.
pub(crate) fn validate_point_inverse_uniformization_roundtrip_impl(
    curve: &AnalyticWeierstrassCurve,
    point: &AnalyticCurvePoint,
    period_config: PeriodRecoveryConfig,
    abel_config: AbelJacobiConfig,
    validation_config: PointRoundTripValidationConfig,
) -> Result<PointRoundTripValidationReport, AnalyticCurveError> {
    let period_basis_report = curve.recover_period_basis(period_config)?;

    validate_point_inverse_uniformization_roundtrip_with_periods_impl(
        curve,
        point,
        period_basis_report.periods(),
        abel_config,
        validation_config,
    )
}
