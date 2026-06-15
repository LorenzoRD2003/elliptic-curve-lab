use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, periods::PeriodRecoveryConfig,
    periods::WeierstrassCubicRoots,
};
use crate::fields::complex_approx::ComplexApprox;
use crate::numerics::{
    cube_root_branches, is_near_pure_cubic_regime, primitive_cube_root_of_unity,
};

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CardanoRootRecoveryDiagnostics {
    pub(crate) cardano_product_residual_norm: f64,
    pub(crate) cardano_discriminant: Complex64,
    pub(crate) selected_u_branch_index: usize,
    pub(crate) selected_v_branch_index: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct InternalCubicRootRecovery {
    pub(crate) roots: WeierstrassCubicRoots,
    pub(crate) newton_iterations_used: usize,
    pub(crate) validation_residual_norm: f64,
    pub(crate) cardano_root_recovery_diagnostics: Option<CardanoRootRecoveryDiagnostics>,
}

impl AnalyticWeierstrassCurve {
    pub(crate) fn recover_weierstrass_cubic_roots_internal(
        &self,
        config: PeriodRecoveryConfig,
    ) -> Result<InternalCubicRootRecovery, AnalyticCurveError> {
        let four = Complex64::new(4.0, 0.0);
        let p = -*self.g2() / four;
        let q = -*self.g3() / four;

        let half_q = q / Complex64::new(2.0, 0.0);
        let cardano_discriminant = half_q.powu(2) + (p / Complex64::new(3.0, 0.0)).powu(3);
        let rough_recovery = Self::choose_rough_root_recovery(p, q, cardano_discriminant, config)?;
        let mut polished_roots = Vec::with_capacity(3);
        let mut total_newton_iterations_used = 0;
        for root in rough_recovery.rough_roots {
            let polished = self.polish_root_with_newton(root, config)?;
            total_newton_iterations_used += polished.iterations_used;
            polished_roots.push(polished.root);
        }

        let roots = WeierstrassCubicRoots::new(
            polished_roots[0],
            polished_roots[1],
            polished_roots[2],
            config.tolerance(),
        )?;

        if !ComplexApprox::is_zero_with_tolerance(
            &roots.x_squared_coefficient(),
            config.tolerance(),
        ) {
            return Err(AnalyticCurveError::CubicRootRecoveryFailed);
        }

        if !ComplexApprox::eq_with_tolerance(&roots.g2(), self.g2(), config.tolerance())
            || !ComplexApprox::eq_with_tolerance(&roots.g3(), self.g3(), config.tolerance())
        {
            return Err(AnalyticCurveError::CubicRootRecoveryFailed);
        }

        let validation_residual_norm = roots
            .x_squared_coefficient()
            .norm()
            .max((roots.g2() - *self.g2()).norm())
            .max((roots.g3() - *self.g3()).norm());

        Ok(InternalCubicRootRecovery {
            roots,
            newton_iterations_used: total_newton_iterations_used,
            validation_residual_norm,
            cardano_root_recovery_diagnostics: rough_recovery.cardano_diagnostics.map(
                |diagnostics| CardanoRootRecoveryDiagnostics {
                    cardano_product_residual_norm: diagnostics.product_residual_norm,
                    cardano_discriminant: diagnostics.discriminant,
                    selected_u_branch_index: diagnostics.u_branch_index,
                    selected_v_branch_index: diagnostics.v_branch_index,
                },
            ),
        })
    }

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
        let branch_choice = choose_consistent_cardano_branches(
            &u_candidates,
            &v_candidates,
            target_product,
            config,
        )?;
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

    fn polish_root_with_newton(
        &self,
        initial_root: Complex64,
        config: PeriodRecoveryConfig,
    ) -> Result<NewtonPolishResult, AnalyticCurveError> {
        let mut root = initial_root;

        for iteration in 0..config.newton_max_iterations() {
            let residual = Complex64::new(4.0, 0.0) * root.powu(3) - *self.g2() * root - *self.g3();
            let derivative = Complex64::new(12.0, 0.0) * root.powu(2) - *self.g2();

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
                    Complex64::new(4.0, 0.0) * root.powu(3) - *self.g2() * root - *self.g3();
                if ComplexApprox::is_zero_with_tolerance(&final_residual, config.tolerance()) {
                    return Ok(NewtonPolishResult {
                        root,
                        iterations_used: iteration + 1,
                    });
                }
            }
        }

        let residual = Complex64::new(4.0, 0.0) * root.powu(3) - *self.g2() * root - *self.g3();
        if ComplexApprox::is_zero_with_tolerance(&residual, config.tolerance()) {
            Ok(NewtonPolishResult {
                root,
                iterations_used: config.newton_max_iterations(),
            })
        } else {
            Err(AnalyticCurveError::CubicRootRecoveryFailed)
        }
    }
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
