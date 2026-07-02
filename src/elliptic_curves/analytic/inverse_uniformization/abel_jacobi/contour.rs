use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    inverse_uniformization::abel_jacobi::config::{AbelJacobiConfig, LegendreContourStrategy},
};
use crate::numerics::{ComplexLineSegment, ComplexRay};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct AbelJacobiContourPlan {
    pub(super) segment: ComplexLineSegment,
    pub(super) ray: ComplexRay,
    pub(super) ray_compact_parameter_max: f64,
    pub(super) anchor_radius: f64,
    pub(super) tail_length: f64,
    pub(super) min_distance_to_branch_points: f64,
}

/// Chooses the deterministic `segment + ray` contour used by the current
/// Abel-Jacobi quadrature in the Legendre `X`-plane.
///
/// Mathematical role:
/// once the curve has been reduced to `Y² = X(X - 1)(X - λ)`, the inverse
/// map is approximated by integrating from the input coordinate `X = x`
/// out toward `∞`. The implementation does not search over arbitrary
/// homotopy classes. Instead, it chooses one simple representative path
/// with two pieces:
///
/// 1. a straight segment from the start point `X` to an anchor point
///    `A = R e^{iθ}`
/// 2. an outgoing ray from `A` in the same direction `e^{iθ}`, later
///    compactified for numerical quadrature
///
/// The radius `R` is chosen as a coarse scale dominating the visible
/// branch geometry near the start point and the singular locus `{0, 1, λ}`:
/// `R = 4 * scale + 2`, where `scale` is the maximum of
/// `|X|`, `|λ|`, and `1`.
///
/// The omitted distances `|X - 1|` and `|X - λ|` are not essential here.
/// By the triangle inequality they are already controlled, up to a fixed
/// constant factor, by `|X|`, `|λ|`, and `1`. Since `R` is only a coarse
/// anchor heuristic, the simpler scale captures the same growth regime while
/// staying easier to explain.
///
/// Why this shape:
/// - the straight segment gives a simple deterministic way to leave the
///   neighborhood of the input point
/// - the ray gives a simple improper-tail model, compatible with the
///   compactification `r = s / (1 - s)`
/// - using the same angle for the segment endpoint and the outgoing ray
///   avoids introducing an extra bend parameter
///
/// Angle selection policy:
/// the current code evaluates a fixed shortlist of candidate angles:
/// the coordinate axes, the argument of the start point, and a
/// few small perturbations of that argument.
/// For each candidate, it samples the segment and the compactified ray,
/// computes the minimum sampled distance to the branch locus
/// `{0, 1, λ}`, and keeps the candidate maximizing that clearance.
///
/// The heuristic is to prefer contours that stay as far as possible from
/// branch points and therefore make branch continuation of `sqrt(X(X-1)(X-λ))`
/// less ambiguous.
///
/// Failure mode:
/// if no candidate has a sampled clearance larger than the explicit
/// proximity threshold, the function returns
/// `AnalyticCurveError::BranchChoiceAmbiguous`. This is mathematically
/// honest: the current contour family is then judged too close to the
/// singular locus for reliable branch tracking.
///
/// Complexity:
/// `Θ(s + r)` per candidate angle, where `s = config.segment_samples` and
/// `r = config.ray_samples`. Since the candidate-angle list has fixed
/// size, the overall complexity remains `Θ(s + r)`.
/// Chooses a Legendre-side contour according to the configured path strategy.
pub(super) fn choose_legendre_contour(
    start_x: Complex64,
    lambda: Complex64,
    config: AbelJacobiConfig,
) -> Result<AbelJacobiContourPlan, AnalyticCurveError> {
    match config.legendre_contour_strategy() {
        LegendreContourStrategy::CanonicalSegmentThenRay => {
            choose_canonical_legendre_contour(start_x, lambda, config)
        }
    }
}

fn choose_canonical_legendre_contour(
    start_x: Complex64,
    lambda: Complex64,
    config: AbelJacobiConfig,
) -> Result<AbelJacobiContourPlan, AnalyticCurveError> {
    let scale = start_x.norm().max(lambda.norm()).max(1.0);
    let anchor_radius = 4.0 * scale + 2.0;
    let tail_length = 4.0 * anchor_radius;
    let branch_locus = [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), lambda];
    let mut best_plan: Option<(AbelJacobiContourPlan, f64)> = None;

    for angle in candidate_ray_angles(start_x) {
        if !angle.is_finite() {
            continue;
        }

        let direction = Complex64::new(angle.cos(), angle.sin());
        let anchor = direction * anchor_radius;
        let segment = ComplexLineSegment::new(start_x, anchor);
        let ray = ComplexRay::new(anchor, angle);
        if !segment.is_finite() || !ray.is_finite() {
            continue;
        }
        let ray_compact_parameter_max = ray.compact_parameter_from_distance(tail_length);

        let min_distance = segment
            .sample_uniform(config.segment_samples())
            .into_iter()
            .chain(ray.sample_compact_parameter(ray_compact_parameter_max, config.ray_samples()))
            .map(|sample| {
                branch_locus
                    .iter()
                    .map(|branch_point| (sample - *branch_point).norm())
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(f64::INFINITY, f64::min);

        let plan = AbelJacobiContourPlan {
            segment,
            ray,
            ray_compact_parameter_max,
            anchor_radius,
            tail_length,
            min_distance_to_branch_points: min_distance,
        };

        let is_better = best_plan
            .as_ref()
            .map(|(_, best_distance)| min_distance > *best_distance)
            .unwrap_or(true);

        if is_better {
            best_plan = Some((plan, min_distance));
        }
    }

    let Some((plan, min_distance)) = best_plan else {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    };

    if min_distance <= contour_proximity_threshold(config) {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    }

    Ok(plan)
}

fn candidate_ray_angles(start_x: Complex64) -> Vec<f64> {
    let base_angle = if start_x == Complex64::new(0.0, 0.0) {
        0.0
    } else {
        start_x.arg()
    };

    vec![
        0.0,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        -std::f64::consts::FRAC_PI_2,
        base_angle,
        base_angle + std::f64::consts::FRAC_PI_4,
        base_angle - std::f64::consts::FRAC_PI_4,
        base_angle + std::f64::consts::FRAC_PI_6,
        base_angle - std::f64::consts::FRAC_PI_6,
    ]
}

fn contour_proximity_threshold(config: AbelJacobiConfig) -> f64 {
    10.0 * config
        .tolerance()
        .absolute
        .max(config.tolerance().relative)
        .max(f64::EPSILON.sqrt())
}

#[cfg(test)]
mod tests {

    use num_complex::Complex64;

    use crate::elliptic_curves::analytic::inverse_uniformization::abel_jacobi::{
        AbelJacobiConfig, AnalyticCurveError, LegendreContourStrategy, choose_legendre_contour,
    };
    use crate::numerics::ApproxTolerance;

    #[test]
    fn contour_strategy_can_report_branch_choice_ambiguity_when_clearance_threshold_is_too_large() {
        let config = AbelJacobiConfig::strict()
            .with_tolerance(ApproxTolerance::new(100.0, 100.0))
            .unwrap()
            .with_max_lattice_corrections(4)
            .unwrap()
            .with_legendre_contour_strategy(LegendreContourStrategy::CanonicalSegmentThenRay)
            .unwrap();

        assert!(matches!(
            choose_legendre_contour(Complex64::new(0.2, 0.15), Complex64::new(0.5, 0.0), config),
            Err(AnalyticCurveError::BranchChoiceAmbiguous)
        ));
    }
}
