use crate::elliptic_curves::{
    CanonicalTauRecoveryReport, CubicRootConfigurationReport, CubicRootRecoveryReport,
    LegendreParameter, LegendreParameterConditioning, LegendreParameterOrbit, LegendreReduction,
    LegendreReductionReport, NumericalRecoveryMetadata, PeriodBasisRecoveryReport,
    PeriodLatticeApprox, PeriodRecoveryConfig, PeriodRecoveryReport,
    RecoveredPeriodBasis, RecoveredPeriodBasisReport, TauRecoveryReport, WeierstrassCubicRoots,
    cubic_root_configuration_report,
};
use crate::visualization::Visualizable;

use super::formatting::{
    format_analytic_cubic_model, format_complex_scalar_compact, format_cubic_root_configuration,
    format_cubic_root_separation, format_legendre_orbit_element_kind,
    format_legendre_parameter_conditioning, format_legendre_scalar,
    format_period_recovery_method, format_period_recovery_status, format_root_scalar,
    roots_need_diagnostic_precision,
};

/// Describes the numerical-policy bundle used for period recovery.
pub fn describe_period_recovery_config(config: &PeriodRecoveryConfig) -> String {
    [
        "Period recovery config".to_string(),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            config.tolerance().absolute,
            config.tolerance().relative
        ),
        format!("Newton iteration budget = {}", config.newton_max_iterations()),
        format!("AGM iteration budget = {}", config.agm_max_iterations()),
        format!(
            "Abel-Jacobi integration steps = {}",
            config.abel_jacobi_integration_steps()
        ),
        format!(
            "branch lattice search radius = {}",
            config.branch_lattice_search_radius()
        ),
        format!(
            "fundamental-domain reduction steps = {}",
            config.fundamental_domain_reduction_max_steps()
        ),
    ]
    .join("\n")
}

/// Describes one chosen approximate period basis.
pub fn describe_period_lattice(periods: &PeriodLatticeApprox) -> String {
    [
        "Approximate period lattice".to_string(),
        format!("ω₁ ≈ {}", format_complex_scalar_compact(periods.omega1())),
        format!("ω₂ ≈ {}", format_complex_scalar_compact(periods.omega2())),
        format!(
            "τ = ω₂ / ω₁ ≈ {}",
            format_complex_scalar_compact(periods.tau().tau())
        ),
        "This is one chosen ordered basis, not a canonical lattice representative.".to_string(),
    ]
    .join("\n")
}

/// Describes one numerical period-recovery metadata bundle.
pub fn describe_numerical_recovery_metadata(metadata: &NumericalRecoveryMetadata) -> String {
    let mut lines = vec![
        "Numerical recovery metadata".to_string(),
        format!(
            "resolved method = {}",
            format_period_recovery_method(metadata.resolved_method())
        ),
        format!("status = {}", format_period_recovery_status(metadata.status())),
        format!(
            "newton iterations used = {}",
            metadata.newton_iterations_used()
        ),
        format!("AGM iterations used = {}", metadata.agm_iterations_used()),
        format!(
            "integration steps used = {}",
            metadata.integration_steps_used()
        ),
        format!(
            "branch lattice searches used = {}",
            metadata.branch_lattice_searches_used()
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            metadata.tolerance().absolute,
            metadata.tolerance().relative
        ),
    ];

    if let Some(residual) = metadata.validation_residual_norm() {
        lines.push(format!("validation residual norm = {:.6e}", residual));
    } else {
        lines.push("validation residual norm = unavailable".to_string());
    }

    if let Some(discriminant) = metadata.cardano_discriminant() {
        lines.push(format!(
            "Cardano discriminant ≈ {}",
            format_complex_scalar_compact(discriminant)
        ));
    }

    if let Some(residual) = metadata.cardano_product_residual_norm() {
        lines.push(format!(
            "Cardano branch product residual norm = {:.6e}",
            residual
        ));
    }

    if let (Some(u_index), Some(v_index)) = (
        metadata.selected_u_branch_index(),
        metadata.selected_v_branch_index(),
    ) {
        lines.push(format!(
            "selected Cardano branch indices = (u: {}, v: {})",
            u_index, v_index
        ));
    }

    if let Some(used_principal) = metadata.used_principal_cardano_branches() {
        lines.push(format!(
            "used principal Cardano branches = {}",
            if used_principal { "yes" } else { "no" }
        ));
    }

    lines.join("\n")
}

/// Describes one validated triple of Weierstrass cubic roots.
pub fn describe_weierstrass_cubic_roots(roots: &WeierstrassCubicRoots) -> String {
    let [first, second, third] = roots.roots();
    let use_diagnostic_precision = roots_need_diagnostic_precision(roots);

    [
        "Weierstrass cubic roots".to_string(),
        format!("root[0] ≈ {}", format_root_scalar(first, use_diagnostic_precision)),
        format!("root[1] ≈ {}", format_root_scalar(second, use_diagnostic_precision)),
        format!("root[2] ≈ {}", format_root_scalar(third, use_diagnostic_precision)),
        "stored order is preserved from construction time and is not canonical".to_string(),
        format!(
            "e₁ + e₂ + e₃ ≈ {}",
            format_root_scalar(&roots.sum(), use_diagnostic_precision)
        ),
        format!(
            "e₁e₂ + e₁e₃ + e₂e₃ ≈ {}",
            format_root_scalar(&roots.pairwise_products_sum(), use_diagnostic_precision)
        ),
        format!(
            "e₁e₂e₃ ≈ {}",
            format_root_scalar(&roots.product(), use_diagnostic_precision)
        ),
        format!("g₂ ≈ {}", format_root_scalar(&roots.g2(), use_diagnostic_precision)),
        format!("g₃ ≈ {}", format_root_scalar(&roots.g3(), use_diagnostic_precision)),
        format!("minimum pairwise distance = {:.6e}", roots.min_pairwise_distance()),
    ]
    .join("\n")
}

pub fn describe_legendre_parameter(parameter: &LegendreParameter) -> String {
    [
        "Legendre parameter".to_string(),
        format!("lambda ≈ {}", format_legendre_scalar(parameter.lambda())),
        format!(
            "1 - lambda ≈ {}",
            format_legendre_scalar(&parameter.one_minus_lambda())
        ),
        "This stores one chosen representative of a six-element S3 orbit.".to_string(),
        "Near-singularity diagnostics depend on a separately supplied tolerance.".to_string(),
    ]
    .join("\n")
}

pub fn describe_legendre_parameter_orbit(orbit: &LegendreParameterOrbit) -> String {
    let mut lines = vec!["Legendre parameter orbit".to_string()];

    for element in orbit.elements() {
        lines.push(format!(
            "{} ≈ {}",
            format_legendre_orbit_element_kind(element.kind()),
            format_legendre_scalar(element.lambda())
        ));
    }

    lines.push(
        "These six values represent the same Legendre class up to root permutation.".to_string(),
    );
    lines.join("\n")
}

pub fn describe_legendre_parameter_conditioning(
    conditioning: LegendreParameterConditioning,
) -> String {
    [
        "Legendre parameter conditioning".to_string(),
        format!(
            "conditioning = {}",
            format_legendre_parameter_conditioning(conditioning)
        ),
        format!(
            "near singular locus = {}",
            if conditioning.is_near_singular() { "yes" } else { "no" }
        ),
    ]
    .join("\n")
}

pub fn describe_legendre_reduction(reduction: &LegendreReduction) -> String {
    let selected = reduction.selected_root_triple();
    let use_diagnostic_root_precision = roots_need_diagnostic_precision(reduction.roots());

    [
        "Legendre reduction".to_string(),
        format!(
            "lambda ≈ {}",
            format_legendre_scalar(reduction.parameter().lambda())
        ),
        format!(
            "selected permutation = [{}, {}, {}]",
            reduction.selected_permutation()[0],
            reduction.selected_permutation()[1],
            reduction.selected_permutation()[2]
        ),
        format!(
            "selected root triple ≈ [{}, {}, {}]",
            format_root_scalar(selected[0], use_diagnostic_root_precision),
            format_root_scalar(selected[1], use_diagnostic_root_precision),
            format_root_scalar(selected[2], use_diagnostic_root_precision)
        ),
        "This selected root triple is used by the reduction and is not a canonical root ordering."
            .to_string(),
        format!(
            "x = {} + ({}) X",
            format_complex_scalar_compact(&reduction.x_translation()),
            format_complex_scalar_compact(&reduction.x_scale())
        ),
        format!(
            "rhs scale factor ≈ {}",
            format_complex_scalar_compact(&reduction.legendre_rhs_scale_factor())
        ),
        format!(
            "principal sqrt(x scale) ≈ {}",
            format_complex_scalar_compact(&reduction.principal_sqrt_x_scale())
        ),
        format!(
            "principal y scale ≈ {}",
            format_complex_scalar_compact(&reduction.legendre_y_scale())
        ),
        format!(
            "invariant differential scale ≈ {}",
            format_complex_scalar_compact(&reduction.invariant_differential_scale())
        ),
        "The principal branch is chosen on sqrt(e1 - e2); flipping it changes only a global sign."
            .to_string(),
    ]
    .join("\n")
}

pub fn describe_legendre_reduction_report(report: &LegendreReductionReport) -> String {
    [
        "Legendre reduction report".to_string(),
        format!("lambda ≈ {}", format_legendre_scalar(report.parameter().lambda())),
        format!(
            "selected orbit element relative to input order = {}",
            format_legendre_orbit_element_kind(
                report.selected_orbit_element_relative_to_input_order()
            )
        ),
        format!(
            "conditioning = {}",
            format_legendre_parameter_conditioning(report.conditioning())
        ),
        format!(
            "near singular locus = {}",
            if report.is_near_singular() { "yes" } else { "no" }
        ),
        format!("singularity distance = {:.6e}", report.singularity_distance()),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
        "This representative was preferred by maximizing min(|lambda|, |1 - lambda|, 1 / |lambda|), then applying deterministic tie-breaks.".to_string(),
        "The orbit-element label is relative to the caller-supplied root order, not canonical by itself.".to_string(),
        format!("reduction summary = {}", report.reduction().format_compact()),
    ]
    .join("\n")
}

pub fn describe_cubic_root_configuration_report(report: &CubicRootConfigurationReport) -> String {
    let mut lines = vec![
        "Cubic-root configuration".to_string(),
        format!(
            "configuration = {}",
            format_cubic_root_configuration(report.configuration())
        ),
        format!(
            "separation = {}",
            format_cubic_root_separation(report.separation())
        ),
        format!("minimum pairwise distance = {:.6e}", report.min_pairwise_distance()),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ];

    match report.conjugate_pair_residual() {
        Some(residual) => lines.push(format!("best conjugate-pair residual = {:.6e}", residual)),
        None => lines.push("best conjugate-pair residual = not applicable".to_string()),
    }

    lines.push(format!("roots summary: {}", report.roots().format_compact()));
    lines.join("\n")
}

pub fn describe_period_recovery_report(report: &PeriodRecoveryReport) -> String {
    [
        "Period recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("ω₁ ≈ {}", format_complex_scalar_compact(report.periods().omega1())),
        format!("ω₂ ≈ {}", format_complex_scalar_compact(report.periods().omega2())),
        format!(
            "τ = ω₂ / ω₁ ≈ {}",
            format_complex_scalar_compact(report.periods().tau().tau())
        ),
        format!(
            "recovered j ≈ {}",
            format_complex_scalar_compact(report.recovered_j())
        ),
        format!(
            "curve-side j ≈ {}",
            format_complex_scalar_compact(report.curve_j())
        ),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(report.difference())
        ),
        format!("|difference| = {:.6e}", report.absolute_difference()),
        format!(
            "agrees under tolerance = {}",
            if report.agrees_approximately() { "yes" } else { "no" }
        ),
    ]
    .join("\n")
}

pub fn describe_recovered_period_basis(basis: &RecoveredPeriodBasis) -> String {
    [
        "Recovered period basis".to_string(),
        format!("ω₁ ≈ {}", format_complex_scalar_compact(basis.omega1())),
        format!("ω₂ ≈ {}", format_complex_scalar_compact(basis.omega2())),
        format!(
            "τ = ω₂ / ω₁ ≈ {}",
            format_complex_scalar_compact(basis.tau().tau())
        ),
        format!("oriented area ≈ {:.6e}", basis.oriented_area()),
        format!("covolume ≈ {:.6e}", basis.covolume()),
        "This is one chosen oriented basis, not a canonical SL₂(ℤ)-class representative."
            .to_string(),
    ]
    .join("\n")
}

pub fn describe_recovered_period_basis_report(report: &RecoveredPeriodBasisReport) -> String {
    [
        "Recovered period basis report".to_string(),
        format!("basis summary = {}", report.basis().format_compact()),
        format!(
            "τ = ω₂ / ω₁ ≈ {}",
            format_complex_scalar_compact(report.tau().tau())
        ),
        format!(
            "invariant differential scale ≈ {}",
            format_complex_scalar_compact(&report.invariant_differential_scale())
        ),
        format!("Legendre reduction summary = {}", report.reduction().format_compact()),
        format!(
            "complete elliptic integral summary = K(λ) ≈ {}, K(1 - λ) ≈ {}",
            format_complex_scalar_compact(report.integral_report().k_lambda.value()),
            format_complex_scalar_compact(report.integral_report().k_complementary.value())
        ),
        "These periods come from transporting the normalized Legendre periods back to the original curve."
            .to_string(),
    ]
    .join("\n")
}

pub fn describe_period_basis_recovery_report(report: &PeriodBasisRecoveryReport) -> String {
    let classification = cubic_root_configuration_report(report.roots(), report.metadata().tolerance());

    [
        "Period-basis recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("roots = {}", report.roots().format_compact()),
        format!(
            "configuration = {}",
            format_cubic_root_configuration(classification.configuration())
        ),
        format!(
            "separation = {}",
            format_cubic_root_separation(classification.separation())
        ),
        format!(
            "Legendre reduction summary = {}",
            report.legendre_reduction().format_compact()
        ),
        format!(
            "λ summary = {}",
            report.legendre_reduction().parameter().format_compact()
        ),
        format!("period basis summary = {}", report.periods().format_compact()),
        format!("τ summary = {}", format_complex_scalar_compact(report.tau().tau())),
        format!("metadata summary = {}", report.metadata().format_compact()),
        "The stored root order is implementation-stable but not canonical.".to_string(),
    ]
    .join("\n")
}

pub fn describe_tau_recovery_report(report: &TauRecoveryReport) -> String {
    [
        "Tau recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("τ ≈ {}", format_complex_scalar_compact(report.tau().tau())),
        format!("period basis summary = {}", report.periods().format_compact()),
        format!("Legendre reduction summary = {}", report.legendre_reduction().format_compact()),
        format!("metadata summary = {}", report.metadata().format_compact()),
        "This τ value is recovered through the full period-basis pipeline, not a separate tau-only algorithm."
            .to_string(),
    ]
    .join("\n")
}

pub fn describe_canonical_tau_recovery_report(report: &CanonicalTauRecoveryReport) -> String {
    [
        "Canonical tau recovery report".to_string(),
        format!(
            "original τ ≈ {}",
            format_complex_scalar_compact(report.original_tau().tau())
        ),
        format!(
            "canonical τ ≈ {}",
            format_complex_scalar_compact(report.canonical_tau().tau())
        ),
        format!(
            "accumulated modular matrix γ = [[{}, {}], [{}, {}]]",
            report.accumulated_matrix().a(),
            report.accumulated_matrix().b(),
            report.accumulated_matrix().c(),
            report.accumulated_matrix().d()
        ),
        format!(
            "fundamental-domain status = {}",
            match report.fundamental_domain_reduction().status() {
                crate::elliptic_curves::FundamentalDomainReductionStatus::AlreadyReduced => {
                    "already reduced"
                }
                crate::elliptic_curves::FundamentalDomainReductionStatus::Reduced => "reduced",
                crate::elliptic_curves::FundamentalDomainReductionStatus::StepLimitReached => {
                    "step limit reached"
                }
            }
        ),
        format!("metadata summary = {}", report.metadata().format_compact()),
        "The canonical τ is obtained by applying the accumulated modular matrix to the naturally recovered τ.".to_string(),
    ]
    .join("\n")
}

pub fn describe_cubic_root_recovery_report(report: &CubicRootRecoveryReport) -> String {
    let classification =
        cubic_root_configuration_report(report.roots(), report.metadata().tolerance());

    [
        "Cubic-root recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("roots = {}", report.roots().format_compact()),
        format!(
            "configuration = {}",
            format_cubic_root_configuration(classification.configuration())
        ),
        format!(
            "separation = {}",
            format_cubic_root_separation(classification.separation())
        ),
        format!(
            "reconstructed g₂ ≈ {}",
            format_complex_scalar_compact(report.reconstructed_g2())
        ),
        format!(
            "curve-side g₂ ≈ {}",
            format_complex_scalar_compact(report.curve_g2())
        ),
        format!(
            "reconstructed g₃ ≈ {}",
            format_complex_scalar_compact(report.reconstructed_g3())
        ),
        format!(
            "curve-side g₃ ≈ {}",
            format_complex_scalar_compact(report.curve_g3())
        ),
        format!("metadata summary = {}", report.metadata().format_compact()),
    ]
    .join("\n")
}

impl Visualizable for PeriodRecoveryConfig {
    fn format_compact(&self) -> String {
        format!(
            "tol=({:.1e}, {:.1e}), Newton≤{}, AGM≤{}, FD≤{}",
            self.tolerance().absolute,
            self.tolerance().relative,
            self.newton_max_iterations(),
            self.agm_max_iterations(),
            self.fundamental_domain_reduction_max_steps()
        )
    }

    fn describe(&self) -> String {
        describe_period_recovery_config(self)
    }
}

impl Visualizable for PeriodLatticeApprox {
    fn format_compact(&self) -> String {
        format!(
            "(ω₁, ω₂) ≈ ({}, {})",
            format_complex_scalar_compact(self.omega1()),
            format_complex_scalar_compact(self.omega2())
        )
    }

    fn describe(&self) -> String {
        describe_period_lattice(self)
    }
}

impl Visualizable for NumericalRecoveryMetadata {
    fn format_compact(&self) -> String {
        format!(
            "{}, {}",
            format_period_recovery_method(self.resolved_method()),
            format_period_recovery_status(self.status())
        )
    }

    fn describe(&self) -> String {
        describe_numerical_recovery_metadata(self)
    }
}

impl Visualizable for LegendreParameter {
    fn format_compact(&self) -> String {
        format!("lambda ≈ {}", format_legendre_scalar(self.lambda()))
    }

    fn describe(&self) -> String {
        describe_legendre_parameter(self)
    }
}

impl Visualizable for LegendreParameterOrbit {
    fn format_compact(&self) -> String {
        let values = self.values();
        format!(
            "[{}, {}, {}, {}, {}, {}]",
            format_legendre_scalar(&values[0]),
            format_legendre_scalar(&values[1]),
            format_legendre_scalar(&values[2]),
            format_legendre_scalar(&values[3]),
            format_legendre_scalar(&values[4]),
            format_legendre_scalar(&values[5]),
        )
    }

    fn describe(&self) -> String {
        describe_legendre_parameter_orbit(self)
    }
}

impl Visualizable for LegendreParameterConditioning {
    fn format_compact(&self) -> String {
        format_legendre_parameter_conditioning(*self).to_string()
    }

    fn describe(&self) -> String {
        describe_legendre_parameter_conditioning(*self)
    }
}

impl Visualizable for LegendreReduction {
    fn format_compact(&self) -> String {
        format!(
            "lambda ≈ {}; perm = [{}, {}, {}]",
            format_legendre_scalar(self.parameter().lambda()),
            self.selected_permutation()[0],
            self.selected_permutation()[1],
            self.selected_permutation()[2],
        )
    }

    fn describe(&self) -> String {
        describe_legendre_reduction(self)
    }
}

impl Visualizable for LegendreReductionReport {
    fn format_compact(&self) -> String {
        format!(
            "{}; {}",
            format_legendre_orbit_element_kind(
                self.selected_orbit_element_relative_to_input_order()
            ),
            format_legendre_parameter_conditioning(self.conditioning())
        )
    }

    fn describe(&self) -> String {
        describe_legendre_reduction_report(self)
    }
}

impl Visualizable for WeierstrassCubicRoots {
    fn format_compact(&self) -> String {
        let [first, second, third] = self.roots();
        format!(
            "[{}, {}, {}]",
            format_complex_scalar_compact(first),
            format_complex_scalar_compact(second),
            format_complex_scalar_compact(third)
        )
    }

    fn describe(&self) -> String {
        describe_weierstrass_cubic_roots(self)
    }
}

impl Visualizable for CubicRootConfigurationReport {
    fn format_compact(&self) -> String {
        format!(
            "{}; {}",
            format_cubic_root_configuration(self.configuration()),
            format_cubic_root_separation(self.separation())
        )
    }

    fn describe(&self) -> String {
        describe_cubic_root_configuration_report(self)
    }
}

impl Visualizable for PeriodRecoveryReport {
    fn format_compact(&self) -> String {
        format!("Δj_recovery ≈ {}", format_complex_scalar_compact(self.difference()))
    }

    fn describe(&self) -> String {
        describe_period_recovery_report(self)
    }
}

impl Visualizable for RecoveredPeriodBasis {
    fn format_compact(&self) -> String {
        format!(
            "(ω₁, ω₂) ≈ ({}, {})",
            format_complex_scalar_compact(self.omega1()),
            format_complex_scalar_compact(self.omega2())
        )
    }

    fn describe(&self) -> String {
        describe_recovered_period_basis(self)
    }
}

impl Visualizable for RecoveredPeriodBasisReport {
    fn format_compact(&self) -> String {
        format!(
            "τ ≈ {}; {}",
            format_complex_scalar_compact(self.tau().tau()),
            self.basis().format_compact()
        )
    }

    fn describe(&self) -> String {
        describe_recovered_period_basis_report(self)
    }
}

impl Visualizable for PeriodBasisRecoveryReport {
    fn format_compact(&self) -> String {
        format!(
            "τ ≈ {}; {}",
            format_complex_scalar_compact(self.tau().tau()),
            self.metadata().format_compact()
        )
    }

    fn describe(&self) -> String {
        describe_period_basis_recovery_report(self)
    }
}

impl Visualizable for TauRecoveryReport {
    fn format_compact(&self) -> String {
        format!("τ ≈ {}", format_complex_scalar_compact(self.tau().tau()))
    }

    fn describe(&self) -> String {
        describe_tau_recovery_report(self)
    }
}

impl Visualizable for CanonicalTauRecoveryReport {
    fn format_compact(&self) -> String {
        format!(
            "τ_FD ≈ {}",
            format_complex_scalar_compact(self.canonical_tau().tau())
        )
    }

    fn describe(&self) -> String {
        describe_canonical_tau_recovery_report(self)
    }
}

impl Visualizable for CubicRootRecoveryReport {
    fn format_compact(&self) -> String {
        format!(
            "Δg₂ ≈ {}, Δg₃ ≈ {}",
            format_complex_scalar_compact(self.g2_comparison().difference()),
            format_complex_scalar_compact(self.g3_comparison().difference())
        )
    }

    fn describe(&self) -> String {
        describe_cubic_root_recovery_report(self)
    }
}
