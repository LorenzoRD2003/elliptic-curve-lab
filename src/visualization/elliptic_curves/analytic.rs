use num_complex::Complex64;

use crate::ComplexApprox;
use crate::elliptic_curves::{
    AnalyticCurveMembershipReport, AnalyticDivisionPolynomialComparisonCase,
    AnalyticDivisionPolynomialComparisonStatus, AnalyticEvenDivisionPolynomialReport,
    AnalyticInvariants, AnalyticOddDivisionPolynomialReport, AnalyticTorsionPointApprox,
    AnalyticWeierstrassCurve, CanonicalTauRecoveryReport, ComplexLattice, CubicRootConfiguration,
    CubicRootConfigurationReport, CubicRootRecoveryReport, CubicRootSeparation,
    EisensteinSumApprox, EllipticFunctionApproximation, EvenDivisionPolynomialVanishingBranch,
    FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
    FundamentalDomainReductionStep, FundamentalDomainReductionStepReason, HasPoleDistance,
    JInvariantComparisonReport, LegendreOrbitElementKind, LegendreParameter,
    LegendreParameterConditioning, LegendreParameterOrbit, LegendreReduction,
    LegendreReductionReport, ModularInvarianceReport, ModularMatrix, ModularQParameter,
    NumericalRecoveryMetadata, PeriodBasisRecoveryReport, PeriodLatticeApprox,
    PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryReport, PeriodRecoveryStatus,
    RecoveredPeriodBasis, RecoveredPeriodBasisReport, ShortWeierstrassCurve, TauRecoveryReport,
    TorusToCurveMapResult, TorusToCurveValues, TruncationConvergenceReport, WeierstrassCubicRoots,
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
    WeierstrassPApprox, WeierstrassPDerivativeApprox, cubic_root_configuration_report,
};
use crate::visualization::Visualizable;
use crate::visualization::elliptic_curves::format_point_compact;
use crate::visualization::fields::{format_complex, format_complex_compact};

fn is_small_real(value: f64) -> bool {
    value.abs() <= 1.0e-12
}

fn is_small_complex(value: &Complex64) -> bool {
    value.norm() <= 1.0e-12
}

fn format_complex_scalar_compact(value: &Complex64) -> String {
    format_complex_compact(value)
}

fn format_decimal_diagnostic(value: f64) -> String {
    let absolute_value = value.abs();
    let mut text = if absolute_value >= 1.0e6 || (absolute_value > 0.0 && absolute_value < 1.0e-6) {
        format!("{value:.12e}")
    } else {
        format!("{value:.12}")
    };

    if text.contains('.') {
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
    }

    if text == "-0" { "0".to_string() } else { text }
}

fn format_complex_scalar_diagnostic(value: &Complex64) -> String {
    if is_small_complex(value) {
        return "0".to_string();
    }

    if is_small_real(value.im) {
        return format_decimal_diagnostic(value.re);
    }

    if is_small_real(value.re) {
        return format!("{}i", format_decimal_diagnostic(value.im));
    }

    let real = format_decimal_diagnostic(value.re);
    let imag = format_decimal_diagnostic(value.im.abs());
    let imag_sign = if value.im < 0.0 { '-' } else { '+' };

    format!("{real} {imag_sign} {imag}i")
}

fn legendre_value_needs_diagnostic_precision(value: &Complex64) -> bool {
    let one = Complex64::new(1.0, 0.0);
    let norm = value.norm();
    let distance_to_one = (*value - one).norm();

    norm >= 1.0e6
        || (norm > 0.0 && norm < 1.0e-4)
        || (distance_to_one > 0.0 && distance_to_one < 1.0e-4)
}

fn format_legendre_scalar(value: &Complex64) -> String {
    if legendre_value_needs_diagnostic_precision(value) {
        format_complex_scalar_diagnostic(value)
    } else {
        format_complex_scalar_compact(value)
    }
}

fn roots_need_diagnostic_precision(roots: &WeierstrassCubicRoots) -> bool {
    let [first, second, third] = roots.roots();
    let formatted = [
        format_complex_scalar_compact(first),
        format_complex_scalar_compact(second),
        format_complex_scalar_compact(third),
    ];

    (formatted[0] == formatted[1] && *first != *second)
        || (formatted[0] == formatted[2] && *first != *third)
        || (formatted[1] == formatted[2] && *second != *third)
}

fn format_root_scalar(value: &Complex64, use_diagnostic_precision: bool) -> String {
    if use_diagnostic_precision {
        format_complex_scalar_diagnostic(value)
    } else {
        format_complex_scalar_compact(value)
    }
}

fn append_polynomial_term(output: &mut String, coefficient: Complex64, suffix: &str) {
    if is_small_complex(&coefficient) {
        return;
    }

    if is_small_real(coefficient.im) {
        if coefficient.re < 0.0 {
            output.push_str(&format!(" - {:.6}{}", coefficient.re.abs(), suffix));
        } else {
            output.push_str(&format!(" + {:.6}{}", coefficient.re, suffix));
        }
    } else {
        output.push_str(&format!(
            " + ({}){}",
            format_complex_scalar_compact(&coefficient),
            suffix
        ));
    }
}

/// Formats the analytic cubic model `y² = 4x³ - g₂x - g₃` while suppressing
/// numerically negligible coefficients.
pub fn format_analytic_cubic_model(curve: &AnalyticWeierstrassCurve) -> String {
    let mut equation = "y^2 = 4x^3".to_string();
    append_polynomial_term(&mut equation, -*curve.g2(), "x");
    append_polynomial_term(&mut equation, -*curve.g3(), "");
    equation
}

/// Formats the short-Weierstrass companion of an analytic curve over the
/// approximate complex backend while suppressing numerically negligible terms.
pub fn format_short_weierstrass_over_complex(
    curve: &ShortWeierstrassCurve<ComplexApprox>,
) -> String {
    let mut equation = "y^2 = x^3".to_string();
    append_polynomial_term(&mut equation, *curve.a(), "x");
    append_polynomial_term(&mut equation, *curve.b(), "");
    equation
}

/// Describes a rank-two complex lattice by its basis and associated shape.
pub fn describe_complex_lattice(lattice: &ComplexLattice) -> String {
    let tau_text = lattice
        .tau()
        .map(|tau| format_complex_scalar_compact(tau.tau()))
        .unwrap_or_else(|_| "unavailable".to_string());

    [
        "Complex lattice".to_string(),
        format!("ω₁ = {}", format_complex_scalar_compact(lattice.omega1())),
        format!("ω₂ = {}", format_complex_scalar_compact(lattice.omega2())),
        format!("τ = ω₂ / ω₁ = {}", tau_text),
        format!("oriented area = {:.6}", lattice.oriented_area()),
        format!("covolume = {:.6}", lattice.covolume()),
    ]
    .join("\n")
}

/// Describes one truncated Eisenstein sum approximation.
pub fn describe_eisenstein_sum(sum: &EisensteinSumApprox) -> String {
    [
        "Eisenstein sum".to_string(),
        format!("weight k = {}", sum.weight),
        format!("truncation radius = {}", sum.truncation.radius()),
        format!("terms used = {}", sum.terms_used),
        format!("value ≈ {}", format_complex_scalar_compact(&sum.value)),
    ]
    .join("\n")
}

/// Describes the modular parameter `q = e^{2π i τ}` attached to one
/// upper-half-plane point `τ`.
pub fn describe_q_parameter(q_parameter: &ModularQParameter) -> String {
    [
        "Modular q-parameter".to_string(),
        format!(
            "τ = {}",
            format_complex_scalar_compact(q_parameter.tau().tau())
        ),
        format!(
            "q = e^(2πiτ) ≈ {}",
            format_complex_scalar_compact(q_parameter.q())
        ),
        format!("|q| = {:.6e}", q_parameter.absolute_value()),
        format!(
            "expected from Im(τ): e^(-2π Im(τ)) = {:.6e}",
            (-std::f64::consts::TAU * q_parameter.tau().imaginary_part()).exp()
        ),
        "Because Im(τ) > 0, this always lies inside the open unit disc.".to_string(),
    ]
    .join("\n")
}

fn format_period_recovery_method(method: PeriodRecoveryMethod) -> &'static str {
    match method {
        PeriodRecoveryMethod::AgmViaLegendre => "AGM via Legendre reduction",
        PeriodRecoveryMethod::NumericalPathIntegral => "numerical path integral",
        PeriodRecoveryMethod::Hybrid => "hybrid",
    }
}

fn format_period_recovery_status(status: PeriodRecoveryStatus) -> &'static str {
    match status {
        PeriodRecoveryStatus::Succeeded => "succeeded",
        PeriodRecoveryStatus::HitIterationLimit => "hit iteration limit",
        PeriodRecoveryStatus::BranchChoiceAmbiguous => "branch choice ambiguous",
        PeriodRecoveryStatus::ValidationFailed => "validation failed",
        PeriodRecoveryStatus::Failed => "failed",
    }
}

fn format_cubic_root_configuration(configuration: CubicRootConfiguration) -> &'static str {
    match configuration {
        CubicRootConfiguration::ThreeApproximatelyReal => "three approximately real",
        CubicRootConfiguration::OneApproximatelyRealTwoApproximatelyConjugate => {
            "one approximately real plus an approximately conjugate pair"
        }
        CubicRootConfiguration::GenericComplex => "generic complex",
    }
}

fn format_cubic_root_separation(separation: CubicRootSeparation) -> &'static str {
    match separation {
        CubicRootSeparation::WellSeparated => "well separated",
        CubicRootSeparation::NearlyRepeated => "nearly repeated",
    }
}

fn format_legendre_orbit_element_kind(kind: LegendreOrbitElementKind) -> &'static str {
    match kind {
        LegendreOrbitElementKind::Lambda => "lambda",
        LegendreOrbitElementKind::OneMinusLambda => "1 - lambda",
        LegendreOrbitElementKind::ReciprocalLambda => "1 / lambda",
        LegendreOrbitElementKind::ReciprocalOneMinusLambda => "1 / (1 - lambda)",
        LegendreOrbitElementKind::LambdaMinusOneOverLambda => "(lambda - 1) / lambda",
        LegendreOrbitElementKind::LambdaOverLambdaMinusOne => "lambda / (lambda - 1)",
    }
}

fn format_legendre_parameter_conditioning(
    conditioning: LegendreParameterConditioning,
) -> &'static str {
    match conditioning {
        LegendreParameterConditioning::Generic => "generic",
        LegendreParameterConditioning::NearZero => "near zero",
        LegendreParameterConditioning::NearOne => "near one",
        LegendreParameterConditioning::NearInfinity => "near infinity",
    }
}

/// Describes the numerical-policy bundle used for period recovery.
pub fn describe_period_recovery_config(config: &PeriodRecoveryConfig) -> String {
    [
        "Period recovery config".to_string(),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            config.tolerance().absolute,
            config.tolerance().relative
        ),
        format!(
            "Newton iteration budget = {}",
            config.newton_max_iterations()
        ),
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
        format!(
            "status = {}",
            format_period_recovery_status(metadata.status())
        ),
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
        format!(
            "root[0] ≈ {}",
            format_root_scalar(first, use_diagnostic_precision)
        ),
        format!(
            "root[1] ≈ {}",
            format_root_scalar(second, use_diagnostic_precision)
        ),
        format!(
            "root[2] ≈ {}",
            format_root_scalar(third, use_diagnostic_precision)
        ),
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
        format!(
            "g₂ ≈ {}",
            format_root_scalar(&roots.g2(), use_diagnostic_precision)
        ),
        format!(
            "g₃ ≈ {}",
            format_root_scalar(&roots.g3(), use_diagnostic_precision)
        ),
        format!(
            "minimum pairwise distance = {:.6e}",
            roots.min_pairwise_distance()
        ),
    ]
    .join("\n")
}

/// Describes one chosen Legendre parameter.
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

/// Describes the six classical transforms in one Legendre orbit.
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

/// Describes the coarse conditioning class of one chosen Legendre parameter.
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
            if conditioning.is_near_singular() {
                "yes"
            } else {
                "no"
            }
        ),
    ]
    .join("\n")
}

/// Describes one affine reduction from a Weierstrass cubic to Legendre form.
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

/// Describes one Legendre reduction report together with its input-order-relative label.
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

/// Describes the geometric configuration and separation status of one cubic-root triple.
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
        format!(
            "minimum pairwise distance = {:.6e}",
            report.min_pairwise_distance()
        ),
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

    lines.push(format!(
        "roots summary: {}",
        report.roots().format_compact()
    ));

    lines.join("\n")
}

/// Describes one side-by-side comparison between the Eisenstein-sum and
/// `q`-expansion routes to the modular `j`-invariant.
pub fn describe_j_invariant_comparison(report: &JInvariantComparisonReport) -> String {
    [
        "j-invariant comparison".to_string(),
        format!("τ = {}", format_complex_scalar_compact(report.tau().tau())),
        format!(
            "lattice truncation radius = {}",
            report.lattice_truncation().radius()
        ),
        format!("q-expansion terms = {}", report.q_truncation().terms()),
        format!(
            "j from Eisenstein sums ≈ {}",
            format_complex_scalar_compact(report.eisenstein_j())
        ),
        format!(
            "j from q-expansion ≈ {}",
            format_complex_scalar_compact(report.q_expansion_j())
        ),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(report.difference())
        ),
        format!("|difference| = {:.6e}", report.absolute_difference()),
        format!(
            "agrees under tolerance = {}",
            if report.agrees_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

/// Describes a recovered-period comparison against the curve-side `j`-invariant.
pub fn describe_period_recovery_report(report: &PeriodRecoveryReport) -> String {
    [
        "Period recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!(
            "ω₁ ≈ {}",
            format_complex_scalar_compact(report.periods().omega1())
        ),
        format!(
            "ω₂ ≈ {}",
            format_complex_scalar_compact(report.periods().omega2())
        ),
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
            if report.agrees_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
    ]
    .join("\n")
}

/// Describes one recovered period basis.
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

/// Describes one Legendre-to-period-basis transport report.
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

/// Describes one full curve-level period-basis recovery report.
pub fn describe_period_basis_recovery_report(report: &PeriodBasisRecoveryReport) -> String {
    let classification =
        cubic_root_configuration_report(report.roots(), report.metadata().tolerance());

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
        format!(
            "period basis summary = {}",
            report.periods().format_compact()
        ),
        format!(
            "τ summary = {}",
            format_complex_scalar_compact(report.tau().tau())
        ),
        format!("metadata summary = {}", report.metadata().format_compact()),
        "The stored root order is implementation-stable but not canonical.".to_string(),
    ]
    .join("\n")
}

/// Describes one `τ`-focused recovery report.
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

/// Describes one canonically normalized `τ` recovery report.
pub fn describe_canonical_tau_recovery_report(report: &CanonicalTauRecoveryReport) -> String {
    [
        "Canonical tau recovery report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
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
                FundamentalDomainReductionStatus::AlreadyReduced => "already reduced",
                FundamentalDomainReductionStatus::Reduced => "reduced",
                FundamentalDomainReductionStatus::StepLimitReached => "step limit reached",
            }
        ),
        format!("metadata summary = {}", report.metadata().format_compact()),
        "The canonical τ is obtained by applying the accumulated modular matrix to the naturally recovered τ."
            .to_string(),
    ]
    .join("\n")
}

/// Describes one cubic-root recovery report together with the reconstruction checks.
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
            "Δg₂ ≈ {}",
            format_complex_scalar_compact(report.g2_comparison().difference())
        ),
        format!(
            "reconstructed g₃ ≈ {}",
            format_complex_scalar_compact(report.reconstructed_g3())
        ),
        format!(
            "curve-side g₃ ≈ {}",
            format_complex_scalar_compact(report.curve_g3())
        ),
        format!(
            "Δg₃ ≈ {}",
            format_complex_scalar_compact(report.g3_comparison().difference())
        ),
        format!(
            "reconstruction agrees under tolerance = {}",
            if report.reconstruction_agrees() {
                "yes"
            } else {
                "no"
            }
        ),
        format!("metadata summary = {}", report.metadata().format_compact()),
    ]
    .join("\n")
}

/// Describes one modular matrix together with its action on the standard
/// generators of the modular group.
pub fn describe_modular_matrix(matrix: &ModularMatrix) -> String {
    [
        "Modular matrix".to_string(),
        format!(
            "γ = [[{}, {}], [{}, {}]]",
            matrix.a(),
            matrix.b(),
            matrix.c(),
            matrix.d()
        ),
        format!("determinant = {}", matrix.determinant()),
        "action on τ: γ(τ) = (aτ + b) / (cτ + d)".to_string(),
    ]
    .join("\n")
}

/// Describes one numerical modular-invariance experiment comparing
/// `j(τ)` and `j(γτ)`.
pub fn describe_modular_invariance_report(report: &ModularInvarianceReport) -> String {
    [
        "Modular invariance check".to_string(),
        format!(
            "original τ = {}",
            format_complex_scalar_compact(report.original_tau().tau())
        ),
        format!(
            "transformed τ = {}",
            format_complex_scalar_compact(report.transformed_tau().tau())
        ),
        format!(
            "γ = [[{}, {}], [{}, {}]]",
            report.matrix().a(),
            report.matrix().b(),
            report.matrix().c(),
            report.matrix().d()
        ),
        format!("truncation radius = {}", report.truncation().radius()),
        format!(
            "j(τ) ≈ {}",
            format_complex_scalar_compact(report.original_j())
        ),
        format!(
            "j(γτ) ≈ {}",
            format_complex_scalar_compact(report.transformed_j())
        ),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(report.difference())
        ),
        format!("|difference| = {:.6e}", report.absolute_difference()),
        format!(
            "invariant under tolerance = {}",
            if report.invariant_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

fn format_fundamental_domain_step_reason(
    reason: FundamentalDomainReductionStepReason,
) -> &'static str {
    match reason {
        FundamentalDomainReductionStepReason::RealPartOutsideCenteredStrip => {
            "real part lay outside the centered strip"
        }
        FundamentalDomainReductionStepReason::NormLessThanOne => "norm was less than one",
    }
}

fn format_fundamental_domain_status(status: FundamentalDomainReductionStatus) -> &'static str {
    match status {
        FundamentalDomainReductionStatus::AlreadyReduced => "already reduced",
        FundamentalDomainReductionStatus::Reduced => "reduced",
        FundamentalDomainReductionStatus::StepLimitReached => "step limit reached",
    }
}

/// Describes one actual modular step applied during reduction to the standard
/// fundamental domain.
pub fn describe_fundamental_domain_reduction_step(step: &FundamentalDomainReductionStep) -> String {
    [
        "Fundamental-domain reduction step".to_string(),
        format!(
            "γ_step = [[{}, {}], [{}, {}]]",
            step.applied_matrix().a(),
            step.applied_matrix().b(),
            step.applied_matrix().c(),
            step.applied_matrix().d()
        ),
        format!(
            "before = {}",
            format_complex_scalar_compact(step.before().tau())
        ),
        format!(
            "after = {}",
            format_complex_scalar_compact(step.after().tau())
        ),
        format!(
            "reason = {}",
            format_fundamental_domain_step_reason(step.reason())
        ),
    ]
    .join("\n")
}

/// Describes one reduction report for the standard fundamental domain of
/// `SL_2(ℤ)`.
pub fn describe_fundamental_domain_reduction_report(
    report: &FundamentalDomainReductionReport,
) -> String {
    let mut lines = vec![
        "Fundamental-domain reduction".to_string(),
        format!(
            "original τ = {}",
            format_complex_scalar_compact(report.original_tau().tau())
        ),
        format!(
            "reduced τ = {}",
            format_complex_scalar_compact(report.reduced_tau().tau())
        ),
        format!(
            "accumulated γ = [[{}, {}], [{}, {}]]",
            report.accumulated_matrix().a(),
            report.accumulated_matrix().b(),
            report.accumulated_matrix().c(),
            report.accumulated_matrix().d()
        ),
        format!("steps used = {}", report.steps().len()),
        format!(
            "status = {}",
            format_fundamental_domain_status(report.status())
        ),
    ];

    if let Some(last_step) = report.steps().last() {
        lines.push(format!(
            "last step reason = {}",
            format_fundamental_domain_step_reason(last_step.reason())
        ));
    }

    lines.join("\n")
}

/// Describes a side-by-side comparison between two Eisenstein truncations.
pub fn describe_truncation_convergence(report: &TruncationConvergenceReport) -> String {
    [
        "Truncation comparison".to_string(),
        format!("small radius = {}", report.small().radius()),
        format!("large radius = {}", report.large().radius()),
        format!(
            "small value ≈ {}",
            format_complex_scalar_compact(report.small_value())
        ),
        format!(
            "large value ≈ {}",
            format_complex_scalar_compact(report.large_value())
        ),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(report.difference())
        ),
        format!("|difference| = {:.6e}", report.absolute_difference()),
    ]
    .join("\n")
}

/// Describes the approximate analytic invariants attached to one lattice.
pub fn describe_analytic_invariants(invariants: &AnalyticInvariants) -> String {
    [
        "Analytic invariants".to_string(),
        format!("truncation radius = {}", invariants.truncation.radius()),
        format!("g₂ ≈ {}", format_complex_scalar_compact(&invariants.g2)),
        format!("g₃ ≈ {}", format_complex_scalar_compact(&invariants.g3)),
        format!(
            "Δ ≈ {}",
            format_complex_scalar_compact(&invariants.discriminant)
        ),
        format!(
            "j ≈ {}",
            format_complex_scalar_compact(&invariants.j_invariant)
        ),
    ]
    .join("\n")
}

/// Describes one approximate curve-membership report in the analytic model
/// `y² = 4x³ - g₂x - g₃`.
pub fn describe_analytic_curve_membership(report: &AnalyticCurveMembershipReport) -> String {
    [
        "Analytic curve membership".to_string(),
        format!("point: {}", format_point_compact(report.point())),
        format!("lhs ≈ {}", format_complex_scalar_compact(report.lhs())),
        format!("rhs ≈ {}", format_complex_scalar_compact(report.rhs())),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(report.difference())
        ),
        format!("|difference| = {:.6e}", report.absolute_error()),
        format!(
            "holds under tolerance = {}",
            if report.is_on_curve() { "yes" } else { "no" }
        ),
    ]
    .join("\n")
}

/// Describes one torus torsion point together with its analytic image on the cubic.
pub fn describe_analytic_torsion_point_approx(point: &AnalyticTorsionPointApprox) -> String {
    [
        "Analytic torsion point".to_string(),
        format!(
            "torus index = ({}, {}; {})",
            point.torus_point().index().a(),
            point.torus_point().index().b(),
            point.torus_point().index().n(),
        ),
        format!(
            "reduced coordinate = ({:.6}, {:.6})",
            point.torus_point().coordinate().u(),
            point.torus_point().coordinate().v(),
        ),
        format!(
            "z = {}",
            format_complex_scalar_compact(point.torus_point().z())
        ),
        format!(
            "curve point = {}",
            format_point_compact(point.curve_point())
        ),
        format!(
            "lies on curve under tolerance = {}",
            if point.lies_on_curve() { "yes" } else { "no" }
        ),
    ]
    .join("\n")
}

fn format_division_polynomial_status(
    status: &AnalyticDivisionPolynomialComparisonStatus,
) -> &'static str {
    match status {
        AnalyticDivisionPolynomialComparisonStatus::PoleAtIdentity => "pole at identity",
        AnalyticDivisionPolynomialComparisonStatus::VanishesApproximately => {
            "vanishes approximately"
        }
        AnalyticDivisionPolynomialComparisonStatus::DoesNotVanishApproximately => {
            "does not vanish approximately"
        }
    }
}

fn format_even_branch(branch: &EvenDivisionPolynomialVanishingBranch) -> &'static str {
    match branch {
        EvenDivisionPolynomialVanishingBranch::YApproxZero => "y(P) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::XCriterionApproxZero => "ε_n(x(P)) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::BothBranches => "both y(P) ≈ 0 and ε_n(x(P)) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::NeitherBranch => {
            "neither y(P) nor ε_n(x(P)) is approximately zero"
        }
    }
}

/// Describes one odd-index analytic torsion comparison through `ψ_n(x)`.
pub fn describe_analytic_odd_division_polynomial_report(
    report: &AnalyticOddDivisionPolynomialReport,
) -> String {
    [
        "Analytic torsion vs division polynomial (odd n)".to_string(),
        format!(
            "torus index = ({}, {}; {})",
            report.torsion_point().torus_point().index().a(),
            report.torsion_point().torus_point().index().b(),
            report.torsion_point().torus_point().index().n(),
        ),
        format!(
            "curve point = {}",
            format_point_compact(report.torsion_point().curve_point())
        ),
        format!(
            "x = ℘(z) ≈ {}",
            format_complex_scalar_compact(report.x_value())
        ),
        format!(
            "ψ_n(x) ≈ {}",
            format_complex_scalar_compact(report.psi_n_x())
        ),
        format!("|ψ_n(x)| = {:.6e}", report.absolute_value()),
        format!(
            "status = {}",
            format_division_polynomial_status(report.status())
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

/// Describes one even-index analytic torsion comparison through `ε_n(x)`.
pub fn describe_analytic_even_division_polynomial_report(
    report: &AnalyticEvenDivisionPolynomialReport,
) -> String {
    [
        "Analytic torsion vs division polynomial (even n)".to_string(),
        format!(
            "torus index = ({}, {}; {})",
            report.torsion_point().torus_point().index().a(),
            report.torsion_point().torus_point().index().b(),
            report.torsion_point().torus_point().index().n(),
        ),
        format!(
            "curve point = {}",
            format_point_compact(report.torsion_point().curve_point())
        ),
        format!(
            "x = ℘(z) ≈ {}",
            format_complex_scalar_compact(report.x_value())
        ),
        format!(
            "ε_n(x) ≈ {}",
            format_complex_scalar_compact(report.epsilon_n_x())
        ),
        format!("|ε_n(x)| = {:.6e}", report.absolute_value()),
        format!("branch = {}", format_even_branch(report.branch())),
        format!(
            "status = {}",
            format_division_polynomial_status(report.status())
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

/// Describes one typed analytic torsion comparison against division polynomials.
pub fn describe_analytic_division_polynomial_comparison(
    report: &AnalyticDivisionPolynomialComparisonCase,
) -> String {
    match report {
        AnalyticDivisionPolynomialComparisonCase::Pole {
            torsion_point,
            tolerance,
        } => [
            "Analytic torsion vs division polynomial".to_string(),
            format!(
                "torus index = ({}, {}; {})",
                torsion_point.torus_point().index().a(),
                torsion_point.torus_point().index().b(),
                torsion_point.torus_point().index().n(),
            ),
            "case = pole at identity".to_string(),
            format!(
                "curve point = {}",
                format_point_compact(torsion_point.curve_point())
            ),
            "no finite x = ℘(z) value is available".to_string(),
            format!(
                "tolerance = abs {:.3e}, rel {:.3e}",
                tolerance.absolute, tolerance.relative
            ),
        ]
        .join("\n"),
        AnalyticDivisionPolynomialComparisonCase::Odd(odd_report) => {
            describe_analytic_odd_division_polynomial_report(odd_report)
        }
        AnalyticDivisionPolynomialComparisonCase::Even(even_report) => {
            describe_analytic_even_division_polynomial_report(even_report)
        }
    }
}

fn describe_elliptic_function_approximation<A>(
    name: &str,
    approximation: &A,
    pole_distance: Option<f64>,
    include_z: bool,
) -> String
where
    A: EllipticFunctionApproximation,
{
    let mut lines = vec![name.to_string()];

    if include_z {
        lines.push(format!(
            "z = {}",
            format_complex_scalar_compact(approximation.z())
        ));
    }

    lines.push(format!(
        "truncation radius = {}",
        approximation.truncation().radius()
    ));
    lines.push(format!("terms used = {}", approximation.terms_used()));
    lines.push(format!(
        "value ≈ {}",
        format_complex_scalar_compact(approximation.value())
    ));

    if let Some(distance) = pole_distance {
        lines.push(format!(
            "nearest inspected pole distance = {:.6e}",
            distance
        ));
    }

    lines.join("\n")
}

/// Describes one truncated `℘(z)` approximation.
pub fn describe_weierstrass_p_approx(approximation: &WeierstrassPApprox) -> String {
    describe_elliptic_function_approximation(
        "Weierstrass ℘ approximation",
        approximation,
        Some(approximation.pole_distance()),
        false,
    )
}

/// Describes one truncated `℘′(z)` approximation.
pub fn describe_weierstrass_p_derivative_approx(
    approximation: &WeierstrassPDerivativeApprox,
) -> String {
    describe_elliptic_function_approximation(
        "Weierstrass ℘′ approximation",
        approximation,
        Some(approximation.pole_distance()),
        false,
    )
}

/// Describes one torus-to-curve mapping result.
pub fn describe_torus_to_curve_map(result: &TorusToCurveMapResult) -> String {
    let mut lines = vec![
        "Torus to curve map".to_string(),
        format!("z = {}", format_complex_scalar_compact(result.z())),
        format!("curve = {}", result.curve()),
        format!("point = {}", format_point_compact(result.point())),
    ];

    match result.values() {
        TorusToCurveValues::Pole => {
            lines.push("values = Pole".to_string());
            lines.push(
                "interpretation: z represents a lattice point, so the map lands at infinity"
                    .to_string(),
            );
        }
        TorusToCurveValues::FiniteValues { p, p_prime } => {
            lines.push(format!("℘(z) ≈ {}", format_complex_scalar_compact(p)));
            lines.push(format!(
                "℘′(z) ≈ {}",
                format_complex_scalar_compact(p_prime)
            ));
        }
    }

    lines.push(format!(
        "lies on curve under tolerance = {}",
        if result.lies_on_curve() { "yes" } else { "no" }
    ));

    lines.join("\n")
}

/// Describes one verification report for
/// `℘′(z)^2 = 4℘(z)^3 - g₂℘(z) - g₃`.
pub fn describe_weierstrass_differential_equation(
    report: &WeierstrassDifferentialEquationReport,
) -> String {
    let mut lines = vec![
        "Weierstrass differential equation".to_string(),
        format!("z = {}", format_complex_scalar_compact(report.z())),
    ];

    match report.values() {
        TorusToCurveValues::Pole => {
            lines.push("values = Pole".to_string());
        }
        TorusToCurveValues::FiniteValues { p, p_prime } => {
            lines.push(format!("℘(z) ≈ {}", format_complex_scalar_compact(p)));
            lines.push(format!(
                "℘′(z) ≈ {}",
                format_complex_scalar_compact(p_prime)
            ));
        }
    }

    lines.push(format!(
        "lhs ≈ {}",
        format_complex_scalar_compact(report.lhs())
    ));
    lines.push(format!(
        "rhs ≈ {}",
        format_complex_scalar_compact(report.rhs())
    ));
    lines.push(format!(
        "difference ≈ {}",
        format_complex_scalar_compact(report.difference())
    ));
    lines.push(format!("|difference| = {:.6e}", report.difference().norm()));
    lines.push(format!(
        "status = {}",
        match report.status() {
            WeierstrassDifferentialEquationStatus::HoldsApproximately => "holds approximately",
            WeierstrassDifferentialEquationStatus::FailsApproximately => "fails approximately",
            WeierstrassDifferentialEquationStatus::Pole => "pole",
        }
    ));
    lines.push(format!(
        "tolerance = abs {:.3e}, rel {:.3e}",
        report.tolerance().absolute,
        report.tolerance().relative
    ));

    lines.join("\n")
}

impl Visualizable for ComplexLattice {
    fn format_compact(&self) -> String {
        format!(
            "Λ = ℤ({}) + ℤ({})",
            format_complex(self.omega1()),
            format_complex(self.omega2())
        )
    }

    fn describe(&self) -> String {
        describe_complex_lattice(self)
    }
}

impl Visualizable for EisensteinSumApprox {
    fn format_compact(&self) -> String {
        format!(
            "G_{}(Λ) ≈ {}",
            self.weight,
            format_complex_scalar_compact(&self.value)
        )
    }

    fn describe(&self) -> String {
        describe_eisenstein_sum(self)
    }
}

impl Visualizable for ModularQParameter {
    fn format_compact(&self) -> String {
        format!("q(τ) ≈ {}", format_complex_scalar_compact(self.q()))
    }

    fn describe(&self) -> String {
        describe_q_parameter(self)
    }
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

impl Visualizable for ModularMatrix {
    fn format_compact(&self) -> String {
        format!(
            "[[{}, {}], [{}, {}]]",
            self.a(),
            self.b(),
            self.c(),
            self.d()
        )
    }

    fn describe(&self) -> String {
        describe_modular_matrix(self)
    }
}

impl Visualizable for JInvariantComparisonReport {
    fn format_compact(&self) -> String {
        format!("Δj ≈ {}", format_complex_scalar_compact(self.difference()))
    }

    fn describe(&self) -> String {
        describe_j_invariant_comparison(self)
    }
}

impl Visualizable for PeriodRecoveryReport {
    fn format_compact(&self) -> String {
        format!(
            "Δj_recovery ≈ {}",
            format_complex_scalar_compact(self.difference())
        )
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

impl Visualizable for ModularInvarianceReport {
    fn format_compact(&self) -> String {
        format!(
            "Δ_mod ≈ {}",
            format_complex_scalar_compact(self.difference())
        )
    }

    fn describe(&self) -> String {
        describe_modular_invariance_report(self)
    }
}

impl Visualizable for FundamentalDomainReductionStep {
    fn format_compact(&self) -> String {
        format!(
            "{} -> {}",
            format_complex_scalar_compact(self.before().tau()),
            format_complex_scalar_compact(self.after().tau())
        )
    }

    fn describe(&self) -> String {
        describe_fundamental_domain_reduction_step(self)
    }
}

impl Visualizable for FundamentalDomainReductionReport {
    fn format_compact(&self) -> String {
        format!(
            "{} -> {}",
            format_complex_scalar_compact(self.original_tau().tau()),
            format_complex_scalar_compact(self.reduced_tau().tau())
        )
    }

    fn describe(&self) -> String {
        describe_fundamental_domain_reduction_report(self)
    }
}

impl Visualizable for TruncationConvergenceReport {
    fn format_compact(&self) -> String {
        format!(
            "Δ_trunc ≈ {}",
            format_complex_scalar_compact(self.difference())
        )
    }

    fn describe(&self) -> String {
        describe_truncation_convergence(self)
    }
}

impl Visualizable for AnalyticInvariants {
    fn format_compact(&self) -> String {
        format!(
            "g₂ ≈ {}, g₃ ≈ {}, j ≈ {}",
            format_complex_scalar_compact(&self.g2),
            format_complex_scalar_compact(&self.g3),
            format_complex_scalar_compact(&self.j_invariant)
        )
    }

    fn describe(&self) -> String {
        describe_analytic_invariants(self)
    }
}

impl Visualizable for WeierstrassPApprox {
    fn format_compact(&self) -> String {
        format!("℘(z) ≈ {}", format_complex_scalar_compact(self.value()))
    }

    fn describe(&self) -> String {
        describe_weierstrass_p_approx(self)
    }
}

impl Visualizable for WeierstrassPDerivativeApprox {
    fn format_compact(&self) -> String {
        format!("℘′(z) ≈ {}", format_complex_scalar_compact(self.value()))
    }

    fn describe(&self) -> String {
        describe_weierstrass_p_derivative_approx(self)
    }
}

impl Visualizable for TorusToCurveMapResult {
    fn format_compact(&self) -> String {
        format!(
            "{} ↦ {}",
            format_complex_scalar_compact(self.z()),
            format_point_compact(self.point())
        )
    }

    fn describe(&self) -> String {
        describe_torus_to_curve_map(self)
    }
}

impl Visualizable for AnalyticTorsionPointApprox {
    fn format_compact(&self) -> String {
        format!(
            "({}, {}; {}) ↦ {}",
            self.torus_point().index().a(),
            self.torus_point().index().b(),
            self.torus_point().index().n(),
            format_point_compact(self.curve_point())
        )
    }

    fn describe(&self) -> String {
        describe_analytic_torsion_point_approx(self)
    }
}

impl Visualizable for WeierstrassDifferentialEquationReport {
    fn format_compact(&self) -> String {
        match self.status() {
            WeierstrassDifferentialEquationStatus::HoldsApproximately => {
                "℘′² = 4℘³ - g₂℘ - g₃ (approx)".to_string()
            }
            WeierstrassDifferentialEquationStatus::FailsApproximately => {
                "℘′² ≠ 4℘³ - g₂℘ - g₃ (approx)".to_string()
            }
            WeierstrassDifferentialEquationStatus::Pole => {
                "℘′² = 4℘³ - g₂℘ - g₃ at a pole".to_string()
            }
        }
    }

    fn describe(&self) -> String {
        describe_weierstrass_differential_equation(self)
    }
}

impl Visualizable for AnalyticOddDivisionPolynomialReport {
    fn format_compact(&self) -> String {
        format!("ψ_n(x) ≈ {}", format_complex_scalar_compact(self.psi_n_x()))
    }

    fn describe(&self) -> String {
        describe_analytic_odd_division_polynomial_report(self)
    }
}

impl Visualizable for AnalyticEvenDivisionPolynomialReport {
    fn format_compact(&self) -> String {
        format!(
            "ε_n(x) ≈ {}",
            format_complex_scalar_compact(self.epsilon_n_x())
        )
    }

    fn describe(&self) -> String {
        describe_analytic_even_division_polynomial_report(self)
    }
}

impl Visualizable for AnalyticDivisionPolynomialComparisonCase {
    fn format_compact(&self) -> String {
        match self {
            AnalyticDivisionPolynomialComparisonCase::Pole { .. } => {
                "division polynomial check at a pole".to_string()
            }
            AnalyticDivisionPolynomialComparisonCase::Odd(report) => report.format_compact(),
            AnalyticDivisionPolynomialComparisonCase::Even(report) => report.format_compact(),
        }
    }

    fn describe(&self) -> String {
        describe_analytic_division_polynomial_comparison(self)
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{
        describe_analytic_division_polynomial_comparison,
        describe_analytic_even_division_polynomial_report, describe_analytic_invariants,
        describe_analytic_odd_division_polynomial_report, describe_analytic_torsion_point_approx,
        describe_canonical_tau_recovery_report, describe_complex_lattice,
        describe_cubic_root_configuration_report, describe_cubic_root_recovery_report,
        describe_eisenstein_sum, describe_fundamental_domain_reduction_report,
        describe_fundamental_domain_reduction_step, describe_j_invariant_comparison,
        describe_legendre_parameter, describe_legendre_parameter_conditioning,
        describe_legendre_parameter_orbit, describe_legendre_reduction,
        describe_legendre_reduction_report, describe_modular_invariance_report,
        describe_modular_matrix, describe_numerical_recovery_metadata,
        describe_period_basis_recovery_report, describe_period_lattice,
        describe_period_recovery_config, describe_period_recovery_report, describe_q_parameter,
        describe_recovered_period_basis, describe_recovered_period_basis_report,
        describe_tau_recovery_report, describe_torus_to_curve_map,
        describe_weierstrass_cubic_roots, describe_weierstrass_differential_equation,
        describe_weierstrass_p_approx, format_analytic_cubic_model, format_complex_scalar_compact,
        format_short_weierstrass_over_complex,
    };
    use crate::elliptic_curves::{
        AnalyticCurvePoint, AnalyticDivisionPolynomialComparisonCase, AnalyticWeierstrassCurve,
        ApproxTolerance, ComplexLattice, EllipticFunctionTruncation, LatticeSumTruncation,
        LegendreParameter, LegendreParameterConditioning, LegendreReduction,
        LegendreReductionReport, ModularMatrix, ModularQParameter, NumericalRecoveryMetadata,
        PeriodLatticeApprox, PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryStatus,
        QExpansionTruncation, RecoveredPeriodBasis, UpperHalfPlanePoint, WeierstrassCubicRoots,
        analytic_invariants, compare_analytic_torsion_with_division_polynomial,
        compare_j_from_eisenstein_and_q_expansion,
        compare_primitive_analytic_torsion_with_division_polynomial,
        cubic_root_configuration_report, g4_sum, map_torus_point_to_curve,
        recover_canonical_tau_from_curve, recover_period_basis, recover_tau_from_curve,
        recover_weierstrass_cubic_roots_with_report, reduce_tau_to_standard_fundamental_domain,
        verify_j_modular_invariance, verify_weierstrass_differential_equation, weierstrass_p,
    };
    use crate::visualization::Visualizable;
    use crate::visualization::elliptic_curves::format_point_compact;

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    #[test]
    fn lattice_description_mentions_basis_and_tau() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_complex_lattice(&lattice);

        assert!(text.contains("Complex lattice"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
    }

    #[test]
    fn eisenstein_description_mentions_weight_and_truncation() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let sum = g4_sum(&lattice, LatticeSumTruncation::default_educational()).unwrap();
        let text = describe_eisenstein_sum(&sum);

        assert!(text.contains("Eisenstein sum"));
        assert!(text.contains("weight k = 4"));
        assert!(text.contains("truncation radius"));
        assert!(text.contains("value"));
    }

    #[test]
    fn analytic_invariant_description_mentions_g2_g3_delta_and_j() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let invariants =
            analytic_invariants(&lattice, LatticeSumTruncation::default_educational()).unwrap();
        let text = describe_analytic_invariants(&invariants);

        assert!(text.contains("Analytic invariants"));
        assert!(text.contains("g₂"));
        assert!(text.contains("g₃"));
        assert!(text.contains("Δ"));
        assert!(text.contains("j"));
    }

    #[test]
    fn q_parameter_description_mentions_tau_q_and_open_unit_disc() {
        let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_q_parameter(&q);

        assert!(text.contains("Modular q-parameter"));
        assert!(text.contains("τ ="));
        assert!(text.contains("q = e^(2πiτ)"));
        assert!(text.contains("|q|"));
        assert!(text.contains("open unit disc"));
    }

    #[test]
    fn period_recovery_config_description_mentions_all_budgets() {
        let config = PeriodRecoveryConfig::educational_default();
        let text = describe_period_recovery_config(&config);

        assert!(text.contains("Period recovery config"));
        assert!(text.contains("Newton iteration budget"));
        assert!(text.contains("AGM iteration budget"));
        assert!(text.contains("Abel-Jacobi integration steps"));
        assert!(text.contains("branch lattice search radius"));
        assert!(text.contains("fundamental-domain reduction steps"));
    }

    #[test]
    fn period_lattice_description_mentions_basis_and_tau() {
        let periods = PeriodLatticeApprox::standard_from_tau(UpperHalfPlanePoint::tau_i());
        let text = describe_period_lattice(&periods);

        assert!(text.contains("Approximate period lattice"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
        assert!(text.contains("not a canonical lattice representative"));
    }

    #[test]
    fn numerical_recovery_metadata_description_mentions_method_status_and_counters() {
        let metadata = NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::ValidationFailed,
            7,
            0,
            0,
            2,
            ApproxTolerance::strict(),
            Some(1.0e-9),
        )
        .with_cardano_diagnostics(Complex64::new(3.0, -4.0), 2.5e-14, 0, 2);
        let text = describe_numerical_recovery_metadata(&metadata);

        assert!(text.contains("Numerical recovery metadata"));
        assert!(text.contains("resolved method = hybrid"));
        assert!(text.contains("status = validation failed"));
        assert!(text.contains("newton iterations used = 7"));
        assert!(text.contains("validation residual norm"));
        assert!(text.contains("Cardano discriminant"));
        assert!(text.contains("Cardano branch product residual norm"));
        assert!(text.contains("selected Cardano branch indices"));
        assert!(text.contains("used principal Cardano branches = no"));
    }

    #[test]
    fn legendre_parameter_description_mentions_lambda_and_orbit_caveat() {
        let parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let text = describe_legendre_parameter(&parameter);

        assert!(text.contains("Legendre parameter"));
        assert!(text.contains("lambda"));
        assert!(text.contains("1 - lambda"));
        assert!(text.contains("six-element S3 orbit"));
    }

    #[test]
    fn legendre_parameter_description_avoids_compact_rounding_near_zero() {
        let parameter = LegendreParameter::new(c(-3.333333333e-8, 0.0)).unwrap();
        let text = describe_legendre_parameter(&parameter);

        assert!(text.contains("lambda ≈ -3.333333333"));
        assert!(text.contains("e-8"));
        assert!(!text.contains("lambda ≈ 0\n"));
    }

    #[test]
    fn legendre_parameter_orbit_description_lists_all_classical_transforms() {
        let parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let text = describe_legendre_parameter_orbit(&parameter.orbit());

        assert!(text.contains("Legendre parameter orbit"));
        assert!(text.contains("lambda ≈"));
        assert!(text.contains("1 - lambda"));
        assert!(text.contains("1 / lambda"));
        assert!(text.contains("same Legendre class"));
    }

    #[test]
    fn legendre_conditioning_description_mentions_singularity_verdict() {
        let text =
            describe_legendre_parameter_conditioning(LegendreParameterConditioning::NearInfinity);

        assert!(text.contains("Legendre parameter conditioning"));
        assert!(text.contains("near infinity"));
        assert!(text.contains("near singular locus = yes"));
    }

    #[test]
    fn legendre_reduction_description_mentions_permutation_and_principal_branch_scales() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let reduction = LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let text = describe_legendre_reduction(&reduction);

        assert!(text.contains("Legendre reduction"));
        assert!(text.contains("selected permutation"));
        assert!(text.contains("not a canonical root ordering"));
        assert!(text.contains("rhs scale factor"));
        assert!(text.contains("principal sqrt(x scale)"));
        assert!(text.contains("principal y scale"));
        assert!(text.contains("invariant differential scale"));
    }

    #[test]
    fn legendre_reduction_report_description_mentions_input_order_caveat() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let report =
            LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let text = describe_legendre_reduction_report(&report);

        assert!(text.contains("Legendre reduction report"));
        assert!(text.contains("selected orbit element relative to input order"));
        assert!(text.contains("conditioning = generic"));
        assert!(text.contains("singularity distance"));
        assert!(text.contains("maximizing min(|lambda|, |1 - lambda|, 1 / |lambda|)"));
        assert!(text.contains("not canonical by itself"));
        assert!(text.contains("reduction summary"));
    }

    #[test]
    fn weierstrass_cubic_roots_description_mentions_roots_and_invariants() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_weierstrass_cubic_roots(&roots);

        assert!(text.contains("Weierstrass cubic roots"));
        assert!(text.contains("root[0]"));
        assert!(text.contains("not canonical"));
        assert!(text.contains("g₂"));
        assert!(text.contains("g₃"));
        assert!(text.contains("minimum pairwise distance"));
    }

    #[test]
    fn weierstrass_cubic_roots_description_uses_diagnostic_precision_for_nearly_colliding_roots() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(1.0 + 1.0e-7, 0.0),
            c(-2.0 - 1.0e-7, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_weierstrass_cubic_roots(&roots);

        assert!(text.contains("root[0] ≈ 1"));
        assert!(text.contains("root[1] ≈ 1.0000001"));
        assert!(text.contains("root[2] ≈ -2.0000001"));
    }

    #[test]
    fn cubic_root_configuration_description_mentions_shape_and_separation() {
        let roots = WeierstrassCubicRoots::new(
            c(2.0, 1.0),
            c(-3.0, 0.0),
            c(2.0, -1.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let report = cubic_root_configuration_report(&roots, ApproxTolerance::strict());
        let text = describe_cubic_root_configuration_report(&report);

        assert!(text.contains("Cubic-root configuration"));
        assert!(text.contains("approximately conjugate pair"));
        assert!(text.contains("separation = well separated"));
        assert!(text.contains("best conjugate-pair residual"));
        assert!(text.contains("roots summary"));
    }

    #[test]
    fn j_invariant_comparison_description_mentions_both_routes_and_difference() {
        let report = compare_j_from_eisenstein_and_q_expansion(
            UpperHalfPlanePoint::tau_i(),
            LatticeSumTruncation::default_educational(),
            QExpansionTruncation::new(3).unwrap(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_j_invariant_comparison(&report);

        assert!(text.contains("j-invariant comparison"));
        assert!(text.contains("j from Eisenstein sums"));
        assert!(text.contains("j from q-expansion"));
        assert!(text.contains("|difference|"));
        assert!(text.contains("agrees under tolerance"));
    }

    #[test]
    fn period_recovery_report_description_mentions_periods_and_j_residual() {
        let tau = UpperHalfPlanePoint::tau_i();
        let curve =
            AnalyticWeierstrassCurve::from_tau(&tau, LatticeSumTruncation::new(12).unwrap())
                .unwrap();
        let periods = PeriodLatticeApprox::standard_from_tau(tau);
        let recovered_j = curve.j_invariant().unwrap();
        let report = crate::elliptic_curves::PeriodRecoveryReport::new(
            curve,
            periods,
            recovered_j,
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_period_recovery_report(&report);

        assert!(text.contains("Period recovery report"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("recovered j"));
        assert!(text.contains("curve-side j"));
        assert!(text.contains("|difference|"));
    }

    #[test]
    fn recovered_period_basis_description_mentions_basis_tau_and_covolume() {
        let basis =
            RecoveredPeriodBasis::new(c(2.0, 0.0), c(1.0, 3.0)).expect("valid positive basis");
        let text = describe_recovered_period_basis(&basis);

        assert!(text.contains("Recovered period basis"));
        assert!(text.contains("ω₁"));
        assert!(text.contains("ω₂"));
        assert!(text.contains("τ = ω₂ / ω₁"));
        assert!(text.contains("covolume"));
        assert!(text.contains("not a canonical"));
    }

    #[test]
    fn recovered_period_basis_report_description_mentions_transport_story() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let period_basis_report = recover_period_basis(
            &AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap(),
            PeriodRecoveryConfig::strict(),
        )
        .unwrap();
        let text = describe_recovered_period_basis_report(period_basis_report.basis_report());

        assert!(text.contains("Recovered period basis report"));
        assert!(text.contains("basis summary"));
        assert!(text.contains("invariant differential scale"));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("complete elliptic integral summary"));
        assert!(text.contains("transporting the normalized Legendre periods"));
        assert!(text.contains("τ = ω₂ / ω₁"));
    }

    #[test]
    fn period_basis_recovery_report_description_mentions_curve_roots_tau_and_metadata() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report = recover_period_basis(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_period_basis_recovery_report(&report);

        assert!(text.contains("Period-basis recovery report"));
        assert!(text.contains("curve ="));
        assert!(text.contains("roots ="));
        assert!(text.contains("configuration ="));
        assert!(text.contains("separation ="));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("period basis summary"));
        assert!(text.contains("τ summary"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("not canonical"));
    }

    #[test]
    fn tau_recovery_report_description_mentions_tau_and_no_parallel_pipeline_caveat() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report = recover_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_tau_recovery_report(&report);

        assert!(text.contains("Tau recovery report"));
        assert!(text.contains("τ ≈"));
        assert!(text.contains("period basis summary"));
        assert!(text.contains("Legendre reduction summary"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("not a separate tau-only algorithm"));
    }

    #[test]
    fn canonical_tau_recovery_report_description_mentions_original_tau_canonical_tau_and_matrix() {
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(0.0, 0.0),
            c(0.5, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let curve = AnalyticWeierstrassCurve::new(roots.g2(), roots.g3()).unwrap();
        let report =
            recover_canonical_tau_from_curve(&curve, PeriodRecoveryConfig::strict()).unwrap();
        let text = describe_canonical_tau_recovery_report(&report);

        assert!(text.contains("Canonical tau recovery report"));
        assert!(text.contains("original τ"));
        assert!(text.contains("canonical τ"));
        assert!(text.contains("accumulated modular matrix"));
        assert!(text.contains("fundamental-domain status"));
        assert!(text.contains("metadata summary"));
        assert!(text.contains("applying the accumulated modular matrix"));
    }

    #[test]
    fn cubic_root_recovery_report_description_mentions_reconstruction_and_metadata() {
        let curve = AnalyticWeierstrassCurve::new(c(28.0, 0.0), c(-24.0, 0.0)).unwrap();
        let report =
            recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
                .unwrap();
        let text = describe_cubic_root_recovery_report(&report);

        assert!(text.contains("Cubic-root recovery report"));
        assert!(text.contains("configuration ="));
        assert!(text.contains("separation ="));
        assert!(text.contains("reconstructed g₂"));
        assert!(text.contains("curve-side g₃"));
        assert!(text.contains("metadata summary"));
    }

    #[test]
    fn modular_matrix_description_mentions_entries_and_action() {
        let text = describe_modular_matrix(&ModularMatrix::s());

        assert!(text.contains("Modular matrix"));
        assert!(text.contains("γ = [[0, -1], [1, 0]]"));
        assert!(text.contains("determinant = 1"));
        assert!(text.contains("action on τ"));
    }

    #[test]
    fn modular_invariance_description_mentions_both_taus_and_difference() {
        let report = verify_j_modular_invariance(
            UpperHalfPlanePoint::tau_i(),
            ModularMatrix::s(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_modular_invariance_report(&report);

        assert!(text.contains("Modular invariance check"));
        assert!(text.contains("original τ ="));
        assert!(text.contains("transformed τ ="));
        assert!(text.contains("j(τ)"));
        assert!(text.contains("j(γτ)"));
        assert!(text.contains("|difference|"));
    }

    #[test]
    fn fundamental_domain_descriptions_mention_status_and_reason() {
        let report = reduce_tau_to_standard_fundamental_domain(
            UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap(),
            8,
        )
        .unwrap();
        let report_text = describe_fundamental_domain_reduction_report(&report);
        let step_text = describe_fundamental_domain_reduction_step(&report.steps()[0]);

        assert!(report_text.contains("Fundamental-domain reduction"));
        assert!(report_text.contains("status = reduced"));
        assert!(report_text.contains("steps used ="));
        assert!(step_text.contains("Fundamental-domain reduction step"));
        assert!(step_text.contains("reason = real part lay outside the centered strip"));
    }

    #[test]
    fn weierstrass_p_description_mentions_pole_distance() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let approximation = weierstrass_p(
            &lattice,
            c(0.2, 0.15),
            EllipticFunctionTruncation::default_educational(),
        )
        .unwrap();
        let text = describe_weierstrass_p_approx(&approximation);

        assert!(text.contains("Weierstrass"));
        assert!(text.contains("truncation radius"));
        assert!(text.contains("nearest inspected pole distance"));
    }

    #[test]
    fn torus_to_curve_description_distinguishes_the_pole_case() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let map = map_torus_point_to_curve(
            &lattice,
            c(0.0, 0.0),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let text = describe_torus_to_curve_map(&map);

        assert!(text.contains("Torus to curve map"));
        assert!(text.contains("values = Pole"));
        assert!(text.contains("infinity"));
    }

    #[test]
    fn analytic_torsion_point_description_mentions_index_z_and_curve_point() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let mapped = crate::elliptic_curves::map_torus_torsion_to_curve(
            &lattice,
            3,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_analytic_torsion_point_approx(&mapped[1]);

        assert!(text.contains("Analytic torsion point"));
        assert!(text.contains("torus index ="));
        assert!(text.contains("z ="));
        assert!(text.contains("curve point ="));
    }

    #[test]
    fn analytic_division_polynomial_description_distinguishes_the_pole_case() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_analytic_division_polynomial_comparison(&reports[0]);

        assert!(text.contains("Analytic torsion vs division polynomial"));
        assert!(text.contains("case = pole at identity"));
        assert!(text.contains("no finite x = ℘(z) value is available"));
    }

    #[test]
    fn odd_division_polynomial_description_mentions_psi_n_and_status() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_primitive_analytic_torsion_with_division_polynomial(
            &lattice,
            3,
            LatticeSumTruncation::larger_for_comparison(),
            EllipticFunctionTruncation::new(6).unwrap(),
            ApproxTolerance::new(1.0e-2, 1.0e-2),
        )
        .unwrap();

        let odd_report = match &reports[0] {
            AnalyticDivisionPolynomialComparisonCase::Odd(odd_report) => odd_report,
            other => panic!("expected odd report, got {other:?}"),
        };
        let text = describe_analytic_odd_division_polynomial_report(odd_report);

        assert!(text.contains("odd n"));
        assert!(text.contains("ψ_n(x)"));
        assert!(text.contains("status ="));
    }

    #[test]
    fn even_division_polynomial_description_mentions_branch_and_epsilon_n() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let reports = compare_primitive_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();

        let even_report = match &reports[0] {
            AnalyticDivisionPolynomialComparisonCase::Even(even_report) => even_report,
            other => panic!("expected even report, got {other:?}"),
        };
        let text = describe_analytic_even_division_polynomial_report(even_report);

        assert!(text.contains("even n"));
        assert!(text.contains("ε_n(x)"));
        assert!(text.contains("branch ="));
        assert!(text.contains("neither y(P) nor ε_n(x(P)) is approximately zero"));
    }

    #[test]
    fn differential_equation_description_mentions_lhs_rhs_and_status() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let text = describe_weierstrass_differential_equation(&report);

        assert!(text.contains("Weierstrass differential equation"));
        assert!(text.contains("lhs"));
        assert!(text.contains("rhs"));
        assert!(text.contains("status"));
    }

    #[test]
    fn visualizable_trait_is_hooked_up_for_analytic_reports() {
        let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
        let q = ModularQParameter::from_tau(UpperHalfPlanePoint::tau_i());
        let period_config = PeriodRecoveryConfig::educational_default();
        let period_lattice = PeriodLatticeApprox::standard_from_tau(UpperHalfPlanePoint::tau_i());
        let metadata = NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            PeriodRecoveryStatus::Succeeded,
            3,
            0,
            0,
            0,
            ApproxTolerance::strict(),
            Some(1.0e-12),
        );
        let roots = WeierstrassCubicRoots::new(
            c(1.0, 0.0),
            c(2.0, 0.0),
            c(-3.0, 0.0),
            ApproxTolerance::strict(),
        )
        .unwrap();
        let legendre_parameter = LegendreParameter::new(c(-0.25, 0.0)).unwrap();
        let legendre_orbit = legendre_parameter.orbit();
        let legendre_reduction =
            LegendreReduction::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let legendre_report =
            LegendreReductionReport::from_roots(&roots, ApproxTolerance::strict()).unwrap();
        let root_configuration = cubic_root_configuration_report(&roots, ApproxTolerance::strict());
        let map = map_torus_point_to_curve(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        let report = verify_weierstrass_differential_equation(
            &lattice,
            c(0.2, 0.15),
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();

        assert!(lattice.format_compact().contains("Λ = ℤ"));
        assert!(q.format_compact().contains("q(τ)"));
        assert!(period_config.describe().contains("Period recovery config"));
        assert!(
            period_lattice
                .describe()
                .contains("Approximate period lattice")
        );
        assert!(metadata.describe().contains("Numerical recovery metadata"));
        assert!(legendre_parameter.describe().contains("Legendre parameter"));
        assert!(
            legendre_orbit
                .describe()
                .contains("Legendre parameter orbit")
        );
        assert!(
            LegendreParameterConditioning::NearInfinity
                .describe()
                .contains("Legendre parameter conditioning")
        );
        assert!(legendre_reduction.describe().contains("Legendre reduction"));
        assert!(
            legendre_report
                .describe()
                .contains("Legendre reduction report")
        );
        assert!(roots.describe().contains("Weierstrass cubic roots"));
        assert!(
            root_configuration
                .describe()
                .contains("Cubic-root configuration")
        );
        assert!(q.describe().contains("Modular q-parameter"));
        assert!(ModularMatrix::s().describe().contains("Modular matrix"));
        assert!(map.describe().contains("Torus to curve map"));
        assert!(
            report
                .describe()
                .contains("Weierstrass differential equation")
        );
        let modular_report = verify_j_modular_invariance(
            UpperHalfPlanePoint::tau_i(),
            ModularMatrix::s(),
            LatticeSumTruncation::larger_for_comparison(),
            ApproxTolerance::strict(),
        )
        .unwrap();
        assert!(
            modular_report
                .describe()
                .contains("Modular invariance check")
        );
        let reduction = reduce_tau_to_standard_fundamental_domain(
            UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap(),
            8,
        )
        .unwrap();
        assert!(
            reduction
                .describe()
                .contains("Fundamental-domain reduction")
        );
        let torsion_comparison = compare_analytic_torsion_with_division_polynomial(
            &lattice,
            2,
            LatticeSumTruncation::default_educational(),
            EllipticFunctionTruncation::default_educational(),
            ApproxTolerance::loose(),
        )
        .unwrap();
        assert!(
            torsion_comparison[0]
                .describe()
                .contains("Analytic torsion vs division polynomial")
        );
        let curve = AnalyticWeierstrassCurve::new(c(28.0, 0.0), c(-24.0, 0.0)).unwrap();
        let cubic_root_report =
            recover_weierstrass_cubic_roots_with_report(&curve, PeriodRecoveryConfig::strict())
                .unwrap();
        assert!(
            cubic_root_report
                .describe()
                .contains("Cubic-root recovery report")
        );
        let infinity = AnalyticCurvePoint::infinity();
        assert_eq!(format_point_compact(&infinity), "O");
    }

    #[test]
    fn specialized_curve_formatters_drop_near_zero_terms_and_imaginary_noise() {
        let analytic =
            AnalyticWeierstrassCurve::new(c(188.94472, -1.0e-15), c(1.0e-15, 2.0e-16)).unwrap();
        let short = analytic.as_short_weierstrass();

        assert_eq!(
            format_analytic_cubic_model(&analytic),
            "y^2 = 4x^3 - 188.944720x"
        );
        assert_eq!(
            format_short_weierstrass_over_complex(&short),
            "y^2 = x^3 - 47.236180x"
        );
    }

    #[test]
    fn compact_complex_formatter_drops_tiny_real_noise_next_to_large_imaginary_part() {
        let value = c(5.0e-7, 60690.762066);

        assert_eq!(format_complex_scalar_compact(&value), "60690.762066i");
        assert_eq!(format_complex_scalar_compact(&c(0.0, 0.0)), "0");
    }
}
