use num_complex::Complex64;

use super::super::{AnalyticCurveError, AnalyticWeierstrassCurve};
use super::{
    NumericalRecoveryMetadata, PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryStatus,
    WeierstrassCubicRoots,
};
use crate::elliptic_curves::analytic::ComplexApproxComparison;
use crate::fields::ComplexApprox;

/// Structured report for one successful cubic-root recovery attempt.
#[derive(Clone, Debug, PartialEq)]
pub struct CubicRootRecoveryReport {
    curve: AnalyticWeierstrassCurve,
    roots: WeierstrassCubicRoots,
    g2_comparison: ComplexApproxComparison,
    g3_comparison: ComplexApproxComparison,
    metadata: NumericalRecoveryMetadata,
}

impl CubicRootRecoveryReport {
    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the recovered cubic roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    /// Returns the comparison between reconstructed and original `g₂`.
    pub fn g2_comparison(&self) -> &ComplexApproxComparison {
        &self.g2_comparison
    }

    /// Returns the comparison between reconstructed and original `g₃`.
    pub fn g3_comparison(&self) -> &ComplexApproxComparison {
        &self.g3_comparison
    }

    /// Returns the numerical execution metadata for the recovery run.
    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        &self.metadata
    }

    /// Returns the reconstructed `g₂` derived from the recovered roots.
    pub fn reconstructed_g2(&self) -> &Complex64 {
        self.g2_comparison.left()
    }

    /// Returns the original curve-side `g₂`.
    pub fn curve_g2(&self) -> &Complex64 {
        self.g2_comparison.right()
    }

    /// Returns the reconstructed `g₃` derived from the recovered roots.
    pub fn reconstructed_g3(&self) -> &Complex64 {
        self.g3_comparison.left()
    }

    /// Returns the original curve-side `g₃`.
    pub fn curve_g3(&self) -> &Complex64 {
        self.g3_comparison.right()
    }

    /// Returns whether both reconstructed coefficients agree approximately
    /// with the original curve-side coefficients.
    pub fn reconstruction_agrees(&self) -> bool {
        self.g2_comparison.agrees_approximately() && self.g3_comparison.agrees_approximately()
    }
}

#[derive(Clone, Debug, PartialEq)]
struct InternalCubicRootRecovery {
    roots: WeierstrassCubicRoots,
    newton_iterations_used: usize,
    validation_residual_norm: f64,
    cardano_product_residual_norm: Option<f64>,
    cardano_discriminant: Option<Complex64>,
    selected_u_branch_index: Option<usize>,
    selected_v_branch_index: Option<usize>,
}

/// Recovers the three roots of the Weierstrass cubic
/// `4x^3 - g₂x - g₃ = 4(x - e₁)(x - e₂)(x - e₃)` directly from the
/// analytic invariants `g₂` and `g₃`.
///
/// This is the invariant-level implementation of the current cubic-root
/// recovery pipeline. It first reconstructs the validated analytic curve
/// `y² = 4x³ - g₂x - g₃`, then applies the same depressed-cubic hybrid
/// recovery and Newton-polishing routine used by
/// [`recover_weierstrass_cubic_roots`].
///
/// Complexity: `Θ(n)`, where `n = config.newton_max_iterations()`.
pub fn recover_weierstrass_cubic_roots_from_invariants(
    g2: &Complex64,
    g3: &Complex64,
    config: PeriodRecoveryConfig,
) -> Result<WeierstrassCubicRoots, AnalyticCurveError> {
    let curve = AnalyticWeierstrassCurve::new(*g2, *g3)?;
    recover_weierstrass_cubic_roots_from_curve_impl(&curve, config).map(|result| result.roots)
}

/// Recovers the three roots of the Weierstrass cubic
/// `4x^3 - g₂x - g₃ = 4(x - e₁)(x - e₂)(x - e₃)`.
///
/// The recovery proceeds in three stages.
///
/// 1. Divide by `4` to obtain the depressed monic cubic
///    `x^3 + px + q = 0` with
///    `p = -g₂ / 4` and `q = -g₃ / 4`.
/// 2. Choose one algebraic seeding strategy:
///    - in the generic regime, apply the classical Cardano ansatz
///      `x = u + v`, where
///      `u^3 = -q/2 + sqrt((q/2)^2 + (p/3)^3)` and
///      `v^3 = -q/2 - sqrt((q/2)^2 + (p/3)^3)`.
///      Because cube roots in `C` are branch-dependent, the implementation
///      enumerates the three branches of each radical and picks a pair
///      whose product satisfies the consistency condition `uv ≈ -p/3`.
///    - in the near-pure-cubic regime, when `|p|` is numerically tiny
///      relative to the natural `|q|^{2/3}` scale, skip that fragile branch
///      check and instead start from the three cube roots of `-q`, namely the
///      exact roots of the limiting equation `x^3 + q = 0`.
/// 3. Polish the three algebraic seeds with Newton iteration on
///    `f(x) = 4x^3 - g₂x - g₃`.
///
/// The final triple is validated by checking:
///
/// - pairwise distinctness under `config.tolerance()`
/// - the depressed-cubic relation `e₁ + e₂ + e₃ ≈ 0`
/// - reconstruction of the original coefficients via
///   `g₂ = -4(e₁e₂ + e₁e₃ + e₂e₃)` and `g₃ = 4e₁e₂e₃`
///
/// The stored root order is the one produced by the chosen algebraic seeding
/// formula above. It is stable for this implementation but does not claim any
/// canonical geometric meaning.
///
/// Complexity: `Θ(n)`, where `n = config.newton_max_iterations()`.
/// The algebraic seeding stage inspects only a constant number of candidates,
/// while the Newton polishing performs at most `3n` updates.
pub fn recover_weierstrass_cubic_roots(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<WeierstrassCubicRoots, AnalyticCurveError> {
    recover_weierstrass_cubic_roots_from_invariants(curve.g2(), curve.g3(), config)
}

/// Recovers the Weierstrass cubic roots together with a structured
/// reconstruction report.
///
/// This uses the same depressed-cubic hybrid algebraic seeding plus
/// Newton-polishing route as
/// [`recover_weierstrass_cubic_roots`], then records the recovered roots and
/// the reconstructed `g₂`/`g₃` comparisons in one report object.
///
/// The embedded [`NumericalRecoveryMetadata`] currently reports:
///
/// - `resolved_method = Hybrid`, since the computation combines algebraic
///   closed-form seeds with numerical Newton polishing
/// - the actual total Newton iterations used across the three roots
/// - `agm_iterations_used = 0`
/// - `integration_steps_used = 0`
/// - `branch_lattice_searches_used = 0`
///
/// Complexity: `Θ(n)`, where `n = config.newton_max_iterations()`.
pub fn recover_weierstrass_cubic_roots_with_report(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<CubicRootRecoveryReport, AnalyticCurveError> {
    let recovery = recover_weierstrass_cubic_roots_from_curve_impl(curve, config)?;
    let g2_comparison =
        ComplexApproxComparison::new(recovery.roots.g2(), *curve.g2(), config.tolerance());
    let g3_comparison =
        ComplexApproxComparison::new(recovery.roots.g3(), *curve.g3(), config.tolerance());
    let mut metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::Hybrid,
        PeriodRecoveryStatus::Succeeded,
        recovery.newton_iterations_used,
        0,
        0,
        0,
        config.tolerance(),
        Some(recovery.validation_residual_norm),
    );

    if let (
        Some(cardano_discriminant),
        Some(cardano_product_residual_norm),
        Some(selected_u_branch_index),
        Some(selected_v_branch_index),
    ) = (
        recovery.cardano_discriminant,
        recovery.cardano_product_residual_norm,
        recovery.selected_u_branch_index,
        recovery.selected_v_branch_index,
    ) {
        metadata = metadata.with_cardano_diagnostics(
            cardano_discriminant,
            cardano_product_residual_norm,
            selected_u_branch_index,
            selected_v_branch_index,
        );
    }

    Ok(CubicRootRecoveryReport {
        curve: curve.clone(),
        roots: recovery.roots,
        g2_comparison,
        g3_comparison,
        metadata,
    })
}

fn recover_weierstrass_cubic_roots_from_curve_impl(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<InternalCubicRootRecovery, AnalyticCurveError> {
    let four = Complex64::new(4.0, 0.0);
    let p = -*curve.g2() / four;
    let q = -*curve.g3() / four;

    let half_q = q / Complex64::new(2.0, 0.0);
    let cardano_discriminant = half_q.powu(2) + (p / Complex64::new(3.0, 0.0)).powu(3);
    let rough_recovery = choose_rough_root_recovery(p, q, cardano_discriminant, config)?;
    let mut polished_roots = Vec::with_capacity(3);
    let mut total_newton_iterations_used = 0;
    for root in rough_recovery.rough_roots {
        let polished = polish_root_with_newton(curve, root, config)?;
        total_newton_iterations_used += polished.iterations_used;
        polished_roots.push(polished.root);
    }

    let roots = WeierstrassCubicRoots::new(
        polished_roots[0],
        polished_roots[1],
        polished_roots[2],
        config.tolerance(),
    )?;

    if !ComplexApprox::is_zero_with_tolerance(&roots.x_squared_coefficient(), config.tolerance()) {
        return Err(AnalyticCurveError::CubicRootRecoveryFailed);
    }

    if !ComplexApprox::eq_with_tolerance(&roots.g2(), curve.g2(), config.tolerance())
        || !ComplexApprox::eq_with_tolerance(&roots.g3(), curve.g3(), config.tolerance())
    {
        return Err(AnalyticCurveError::CubicRootRecoveryFailed);
    }

    let validation_residual_norm = roots
        .x_squared_coefficient()
        .norm()
        .max((roots.g2() - *curve.g2()).norm())
        .max((roots.g3() - *curve.g3()).norm());

    Ok(InternalCubicRootRecovery {
        roots,
        newton_iterations_used: total_newton_iterations_used,
        validation_residual_norm,
        cardano_product_residual_norm: rough_recovery
            .cardano_diagnostics
            .map(|diagnostics| diagnostics.product_residual_norm),
        cardano_discriminant: rough_recovery
            .cardano_diagnostics
            .map(|diagnostics| diagnostics.discriminant),
        selected_u_branch_index: rough_recovery
            .cardano_diagnostics
            .map(|diagnostics| diagnostics.u_branch_index),
        selected_v_branch_index: rough_recovery
            .cardano_diagnostics
            .map(|diagnostics| diagnostics.v_branch_index),
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct RoughRootRecovery {
    rough_roots: [Complex64; 3],
    cardano_diagnostics: Option<CardanoDiagnostics>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct CardanoDiagnostics {
    discriminant: Complex64,
    product_residual_norm: f64,
    u_branch_index: usize,
    v_branch_index: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct CardanoBranchChoice {
    u: Complex64,
    v: Complex64,
    u_branch_index: usize,
    v_branch_index: usize,
    product_residual_norm: f64,
}

/// Chooses the algebraic seeds that will later be Newton-polished.
///
/// For generic depressed cubics we use the classical Cardano branch search.
/// When `|p|` is tiny compared with the natural `|q|^{2/3}` scale, the
/// limiting equation is numerically much closer to `x^3 + q = 0` than to a
/// well-conditioned Cardano branch-matching problem, so we seed Newton from
/// the three cube roots of `-q` instead.
fn choose_rough_root_recovery(
    p: Complex64,
    q: Complex64,
    cardano_discriminant: Complex64,
    config: PeriodRecoveryConfig,
) -> Result<RoughRootRecovery, AnalyticCurveError> {
    if is_near_pure_cubic_regime(p, q, config.tolerance()) {
        return Ok(RoughRootRecovery {
            rough_roots: cube_root_branches(-q),
            cardano_diagnostics: None,
        });
    }

    let half_q = q / Complex64::new(2.0, 0.0);
    let radical_sqrt = cardano_discriminant.sqrt();
    let radical_a = -half_q + radical_sqrt;
    let radical_b = -half_q - radical_sqrt;

    let u_candidates = cube_root_branches(radical_a);
    let v_candidates = cube_root_branches(radical_b);
    let target_product = -p / Complex64::new(3.0, 0.0);
    let branch_choice =
        choose_consistent_cardano_branches(&u_candidates, &v_candidates, target_product, config)?;
    let rough_roots = cardano_roots_from_branch_choice(branch_choice.u, branch_choice.v);

    Ok(RoughRootRecovery {
        rough_roots,
        cardano_diagnostics: Some(CardanoDiagnostics {
            discriminant: cardano_discriminant,
            product_residual_norm: branch_choice.product_residual_norm,
            u_branch_index: branch_choice.u_branch_index,
            v_branch_index: branch_choice.v_branch_index,
        }),
    })
}

/// Heuristically detects the regime where the depressed cubic
/// `x^3 + px + q = 0` is numerically close to the pure cubic `x^3 + q = 0`.
///
/// The ratio `|p| / |q|^{2/3}` is the natural dimensionless measure here:
/// if `x` has the pure-cubic scale `|x| ≈ |q|^{1/3}`, then the perturbation
/// term `px` is smaller than the main terms by roughly that factor.
///
/// The tolerance enters through a deliberately coarse relaxation. This is
/// intentionally much looser than the final equality tolerance because the
/// goal is not to certify that `p = 0`, but to detect when strict Cardano
/// branch matching is more fragile than simply seeding Newton from the
/// pure-cubic limit.
fn is_near_pure_cubic_regime(
    p: Complex64,
    q: Complex64,
    tolerance: crate::numerics::ApproxTolerance,
) -> bool {
    let tolerance_scale = tolerance.absolute.max(tolerance.relative).max(f64::EPSILON);
    let q_scale_squared = q.norm().cbrt().powi(2).max(tolerance_scale);
    let normalized_p = p.norm() / q_scale_squared;
    normalized_p <= tolerance_scale.powf(0.25)
}

fn cardano_roots_from_branch_choice(u: Complex64, v: Complex64) -> [Complex64; 3] {
    let omega = primitive_cube_root_of_unity();
    let omega_sq = omega * omega;
    [u + v, omega * u + omega_sq * v, omega_sq * u + omega * v]
}

fn choose_consistent_cardano_branches(
    u_candidates: &[Complex64; 3],
    v_candidates: &[Complex64; 3],
    target_product: Complex64,
    config: PeriodRecoveryConfig,
) -> Result<CardanoBranchChoice, AnalyticCurveError> {
    let mut best_pair: Option<CardanoBranchChoice> = None;

    for (u_branch_index, u) in u_candidates.iter().enumerate() {
        for (v_branch_index, v) in v_candidates.iter().enumerate() {
            let product_residual_norm = (*u * *v - target_product).norm();
            let candidate = CardanoBranchChoice {
                u: *u,
                v: *v,
                u_branch_index,
                v_branch_index,
                product_residual_norm,
            };

            match best_pair {
                None => best_pair = Some(candidate),
                Some(best) if product_residual_norm < best.product_residual_norm => {
                    best_pair = Some(candidate)
                }
                _ => {}
            }
        }
    }

    let Some(best_pair) = best_pair else {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    };

    if !ComplexApprox::eq_with_tolerance(
        &(best_pair.u * best_pair.v),
        &target_product,
        config.tolerance(),
    ) {
        return Err(AnalyticCurveError::BranchChoiceAmbiguous);
    }

    Ok(best_pair)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct NewtonPolishResult {
    root: Complex64,
    iterations_used: usize,
}

fn polish_root_with_newton(
    curve: &AnalyticWeierstrassCurve,
    initial_root: Complex64,
    config: PeriodRecoveryConfig,
) -> Result<NewtonPolishResult, AnalyticCurveError> {
    let mut root = initial_root;

    for iteration in 0..config.newton_max_iterations() {
        let residual = Complex64::new(4.0, 0.0) * root.powu(3) - *curve.g2() * root - *curve.g3();
        let derivative = Complex64::new(12.0, 0.0) * root.powu(2) - *curve.g2();

        if ComplexApprox::is_zero_with_tolerance(&residual, config.tolerance()) {
            return Ok(NewtonPolishResult {
                root,
                iterations_used: iteration,
            });
        }

        if ComplexApprox::is_zero_with_tolerance(&derivative, config.tolerance()) {
            return Err(AnalyticCurveError::CubicRootRecoveryFailed);
        }

        let step = residual / derivative;
        root -= step;

        if ComplexApprox::is_zero_with_tolerance(&step, config.tolerance()) {
            let final_residual =
                Complex64::new(4.0, 0.0) * root.powu(3) - *curve.g2() * root - *curve.g3();
            if ComplexApprox::is_zero_with_tolerance(&final_residual, config.tolerance()) {
                return Ok(NewtonPolishResult {
                    root,
                    iterations_used: iteration + 1,
                });
            }
        }
    }

    let residual = Complex64::new(4.0, 0.0) * root.powu(3) - *curve.g2() * root - *curve.g3();
    if ComplexApprox::is_zero_with_tolerance(&residual, config.tolerance()) {
        Ok(NewtonPolishResult {
            root,
            iterations_used: config.newton_max_iterations(),
        })
    } else {
        Err(AnalyticCurveError::CubicRootRecoveryFailed)
    }
}

fn cube_root_branches(z: Complex64) -> [Complex64; 3] {
    if ComplexApprox::is_zero_with_tolerance(&z, ComplexApprox::default_tolerance()) {
        return [Complex64::new(0.0, 0.0); 3];
    }

    let principal = principal_cube_root(z);
    let omega = primitive_cube_root_of_unity();
    let omega_sq = omega * omega;

    [principal, omega * principal, omega_sq * principal]
}

fn principal_cube_root(z: Complex64) -> Complex64 {
    let radius = z.norm().cbrt();
    let angle = z.arg() / 3.0;
    Complex64::from_polar(radius, angle)
}

fn primitive_cube_root_of_unity() -> Complex64 {
    Complex64::new(-0.5, f64::sqrt(3.0) * 0.5)
}
