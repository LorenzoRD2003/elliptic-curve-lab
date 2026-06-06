use num_complex::Complex64;

use crate::elliptic_curves::analytic::AnalyticCurveError;
use crate::elliptic_curves::analytic::inverse_uniformization::abel_jacobi::config::AbelJacobiConfig;
use crate::elliptic_curves::analytic::inverse_uniformization::abel_jacobi::report::AbelJacobiInitialBranchChoice;
use crate::elliptic_curves::analytic::periods::PeriodRecoveryConfig;
use crate::numerics::{
    ComplexLineSegment, ComplexRay, SimpsonQuadratureDomain,
    composite_simpson_integrate_complex_in_domain,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct BranchState {
    pub(super) sqrt_value: Complex64,
    uses_principal_branch: bool,
}

impl BranchState {
    pub(super) fn initial_branch_choice(self) -> AbelJacobiInitialBranchChoice {
        if self.uses_principal_branch {
            AbelJacobiInitialBranchChoice::Principal
        } else {
            AbelJacobiInitialBranchChoice::Alternate
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ParameterInterval {
    start: f64,
    end: f64,
    subintervals: usize,
}

pub(super) fn derived_root_recovery_config(config: AbelJacobiConfig) -> PeriodRecoveryConfig {
    PeriodRecoveryConfig::new(
        config.tolerance,
        (config.max_branch_adjustments * 2).max(8),
        6,
        config.integration_steps.max(8),
        config.max_lattice_corrections.max(1),
        8,
    )
    .expect("derived Abel-Jacobi root-recovery config must stay valid")
}

/// Chooses the initial sign of `sqrt(X(X-1)(X-λ))` at the starting Legendre
/// point before branch continuation along the contour begins.
///
/// After reducing the curve to `Y² = X(X-1)(X-λ)`, the Abel-Jacobi convention
/// used here is
/// `z = ∫_x^∞ dt / sqrt(4t³ - g₂ t - g₃)`.
/// Under the current change of variables, that means the starting square-root
/// value should be compatible with `-Y`, not with `+Y`.
///
/// Concretely, the code compares the principal square root and its negative
/// against the target value `-transformed_y`, then picks the closer one. If
/// both signs are numerically indistinguishable at the requested tolerance in
/// a regime where the target value is itself too small to break the tie
/// robustly, the function reports `AnalyticCurveError::BranchChoiceAmbiguous`
/// instead of pretending the initial branch was determined reliably.
///
/// Complexity: `Θ(1)`.
pub(super) fn initialize_legendre_branch_state(
    transformed_x: Complex64,
    transformed_y: Complex64,
    lambda: Complex64,
    config: AbelJacobiConfig,
) -> Result<BranchState, AnalyticCurveError> {
    let target_branch_value = -transformed_y;
    let principal = legendre_cubic(transformed_x, lambda).sqrt();
    let alternate = -principal;
    let principal_distance = (principal - target_branch_value).norm();
    let alternate_distance = (alternate - target_branch_value).norm();

    if config
        .tolerance
        .real_close(principal_distance, alternate_distance)
        && target_branch_value.norm()
            <= config
                .tolerance
                .absolute
                .max(config.tolerance.relative * target_branch_value.norm().max(1.0))
    {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    }

    if principal_distance <= alternate_distance {
        Ok(BranchState {
            sqrt_value: principal,
            uses_principal_branch: true,
        })
    } else {
        Ok(BranchState {
            sqrt_value: alternate,
            uses_principal_branch: false,
        })
    }
}

pub(super) fn integrate_segment_with_branch_tracking(
    segment: &ComplexLineSegment,
    lambda: Complex64,
    initial_branch: BranchState,
    subintervals: usize,
    config: AbelJacobiConfig,
) -> Result<(Complex64, BranchState, usize), AnalyticCurveError> {
    integrate_parameterized_path_with_branch_tracking(
        ParameterInterval {
            start: 0.0,
            end: 1.0,
            subintervals,
        },
        |t| segment.point_at(t),
        |_| segment.displacement(),
        lambda,
        initial_branch,
        config,
    )
}

pub(super) fn integrate_ray_with_branch_tracking(
    ray: &ComplexRay,
    s_max: f64,
    lambda: Complex64,
    initial_branch: BranchState,
    subintervals: usize,
    config: AbelJacobiConfig,
) -> Result<(Complex64, BranchState, usize), AnalyticCurveError> {
    integrate_parameterized_path_with_branch_tracking(
        ParameterInterval {
            start: 0.0,
            end: s_max,
            subintervals,
        },
        |s| ray.point_at_compact_parameter(s),
        |s| ray.compact_parameter_derivative(s),
        lambda,
        initial_branch,
        config,
    )
}

/// Integrates one parameterized contour while transporting a square-root
/// branch by continuity from sample to sample.
///
/// The closure pair `point_at`, `derivative_at` defines a complex path
/// `X(u)` and its derivative `dX/du`. At each Simpson node we evaluate
/// `dX/du / √(X(u)(X(u)-1)(X(u)-λ))`, choosing the sign of the square root
/// by minimizing the jump from the previous sample.
///
/// Complexity: `Θ(m)`, where `m` is the normalized Simpson node count for the
/// supplied parameter interval.
fn integrate_parameterized_path_with_branch_tracking<F, G>(
    interval: ParameterInterval,
    point_at: F,
    derivative_at: G,
    lambda: Complex64,
    initial_branch: BranchState,
    config: AbelJacobiConfig,
) -> Result<(Complex64, BranchState, usize), AnalyticCurveError>
where
    F: Fn(f64) -> Complex64,
    G: Fn(f64) -> Complex64,
{
    let domain = SimpsonQuadratureDomain::new(interval.start, interval.end, interval.subintervals)
        .map_err(|_| AnalyticCurveError::AbelJacobiIntegrationFailed)?;
    let mut current_branch = initial_branch;
    let mut branch_adjustments_used = 0usize;

    let integral = composite_simpson_integrate_complex_in_domain(&domain, |index, parameter| {
        let point = point_at(parameter);
        let derivative = derivative_at(parameter);

        if !point.is_finite() || !derivative.is_finite() {
            return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
        }

        if index > 0 {
            let (next_branch, adjusted_here) =
                continue_legendre_branch(current_branch, point, lambda, config)?;
            current_branch = next_branch;
            branch_adjustments_used += usize::from(adjusted_here);

            if branch_adjustments_used > config.max_branch_adjustments {
                return Err(AnalyticCurveError::BranchChoiceAmbiguous);
            }
        }

        Ok(derivative / current_branch.sqrt_value)
    })
    .map_err(AnalyticCurveError::from)?;

    Ok((integral, current_branch, branch_adjustments_used))
}

fn continue_legendre_branch(
    previous_branch: BranchState,
    point: Complex64,
    lambda: Complex64,
    config: AbelJacobiConfig,
) -> Result<(BranchState, bool), AnalyticCurveError> {
    let principal = legendre_cubic(point, lambda).sqrt();
    let alternate = -principal;
    let principal_distance = (principal - previous_branch.sqrt_value).norm();
    let alternate_distance = (alternate - previous_branch.sqrt_value).norm();

    if config
        .tolerance
        .real_close(principal_distance, alternate_distance)
    {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    }

    if principal_distance < alternate_distance {
        let branch = BranchState {
            sqrt_value: principal,
            uses_principal_branch: true,
        };
        Ok((
            branch,
            branch.uses_principal_branch != previous_branch.uses_principal_branch,
        ))
    } else {
        let branch = BranchState {
            sqrt_value: alternate,
            uses_principal_branch: false,
        };
        Ok((
            branch,
            branch.uses_principal_branch != previous_branch.uses_principal_branch,
        ))
    }
}

/// Approximates the omitted improper tail of the Legendre-side contour
/// integral after the sampled ray has been truncated at a finite endpoint.
///
/// If the ray endpoint is `X_end` and the currently continued branch value is
/// `sqrt(X_end(X_end-1)(X_end-λ))`, the implementation uses the first-order
/// asymptotic model
/// `∫_{X_end}^{∞} dX / sqrt(X(X-1)(X-λ)) ≈ 2 X_end / sqrt(X_end(X_end-1)(X_end-λ))`.
///
/// This is the educational leading term coming from the fact that for large
/// `X` the cubic behaves like `X³`, so the integrand behaves like `X^{-3/2}`
/// and its antiderivative behaves like `-2 X^{-1/2}`.
///
/// Complexity: `Θ(1)`.
pub(super) fn estimate_legendre_tail_correction(
    tail_endpoint: Complex64,
    end_branch: BranchState,
) -> Result<Complex64, AnalyticCurveError> {
    if !tail_endpoint.is_finite()
        || !end_branch.sqrt_value.is_finite()
        || tail_endpoint == Complex64::new(0.0, 0.0)
    {
        return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
    }

    Ok(Complex64::new(2.0, 0.0) * tail_endpoint / end_branch.sqrt_value)
}

/// Evaluates the Legendre cubic polynomial `X(X-1)(X-λ)`.
///
/// This is the branch polynomial appearing in the reduced model
/// `Y² = X(X-1)(X-λ)`. The Abel-Jacobi integrand uses its square root, so this
/// helper is the single local source for that cubic expression inside the
/// branch-tracking code.
///
/// Complexity: `Θ(1)`.
fn legendre_cubic(x: Complex64, lambda: Complex64) -> Complex64 {
    x * (x - Complex64::new(1.0, 0.0)) * (x - lambda)
}
