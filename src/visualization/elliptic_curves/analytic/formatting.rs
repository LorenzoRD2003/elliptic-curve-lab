use num_complex::Complex64;

use crate::elliptic_curves::analytic::torsion::{
    AnalyticDivisionPolynomialComparisonStatus, EvenDivisionPolynomialVanishingBranch,
};
use crate::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve, CubicRootConfiguration, CubicRootSeparation,
    FundamentalDomainReductionStatus, FundamentalDomainReductionStepReason,
    InvariantRecoveryInterpretation, LegendreOrbitElementKind, LegendreParameterConditioning,
    PeriodRecoveryMethod, PeriodRecoveryStatus, WeierstrassCubicRoots,
};
use crate::visualization::fields::complex_approx::format_complex_compact;

pub(crate) fn is_small_real(value: f64) -> bool {
    value.abs() <= 1.0e-12
}

pub(crate) fn is_small_complex(value: &Complex64) -> bool {
    value.norm() <= 1.0e-12
}

pub(crate) fn format_complex_scalar_compact(value: &Complex64) -> String {
    format_complex_compact(value)
}

pub(crate) fn format_decimal_diagnostic(value: f64) -> String {
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

pub(crate) fn format_complex_scalar_diagnostic(value: &Complex64) -> String {
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

pub(crate) fn legendre_value_needs_diagnostic_precision(value: &Complex64) -> bool {
    let one = Complex64::new(1.0, 0.0);
    let norm = value.norm();
    let distance_to_one = (*value - one).norm();

    norm >= 1.0e6
        || (norm > 0.0 && norm < 1.0e-4)
        || (distance_to_one > 0.0 && distance_to_one < 1.0e-4)
}

pub(crate) fn format_legendre_scalar(value: &Complex64) -> String {
    if legendre_value_needs_diagnostic_precision(value) {
        format_complex_scalar_diagnostic(value)
    } else {
        format_complex_scalar_compact(value)
    }
}

pub(crate) fn roots_need_diagnostic_precision(roots: &WeierstrassCubicRoots) -> bool {
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

pub(crate) fn format_root_scalar(value: &Complex64, use_diagnostic_precision: bool) -> String {
    if use_diagnostic_precision {
        format_complex_scalar_diagnostic(value)
    } else {
        format_complex_scalar_compact(value)
    }
}

pub(crate) fn append_polynomial_term(output: &mut String, coefficient: Complex64, suffix: &str) {
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

pub(crate) fn format_period_recovery_method(method: PeriodRecoveryMethod) -> &'static str {
    match method {
        PeriodRecoveryMethod::AgmViaLegendre => "AGM via Legendre reduction",
        PeriodRecoveryMethod::NumericalPathIntegral => "numerical path integral",
        PeriodRecoveryMethod::Hybrid => "hybrid",
    }
}

pub(crate) fn format_period_recovery_status(status: PeriodRecoveryStatus) -> &'static str {
    match status {
        PeriodRecoveryStatus::Succeeded => "succeeded",
        PeriodRecoveryStatus::HitIterationLimit => "hit iteration limit",
        PeriodRecoveryStatus::BranchChoiceAmbiguous => "branch choice ambiguous",
        PeriodRecoveryStatus::ValidationFailed => "validation failed",
        PeriodRecoveryStatus::Failed => "failed",
    }
}

pub(crate) fn format_cubic_root_configuration(
    configuration: CubicRootConfiguration,
) -> &'static str {
    match configuration {
        CubicRootConfiguration::ThreeApproximatelyReal => "three approximately real",
        CubicRootConfiguration::OneApproximatelyRealTwoApproximatelyConjugate => {
            "one approximately real plus an approximately conjugate pair"
        }
        CubicRootConfiguration::GenericComplex => "generic complex",
    }
}

pub(crate) fn format_cubic_root_separation(separation: CubicRootSeparation) -> &'static str {
    match separation {
        CubicRootSeparation::WellSeparated => "well separated",
        CubicRootSeparation::NearlyRepeated => "nearly repeated",
    }
}

pub(crate) fn format_legendre_orbit_element_kind(kind: LegendreOrbitElementKind) -> &'static str {
    match kind {
        LegendreOrbitElementKind::Lambda => "lambda",
        LegendreOrbitElementKind::OneMinusLambda => "1 - lambda",
        LegendreOrbitElementKind::ReciprocalLambda => "1 / lambda",
        LegendreOrbitElementKind::ReciprocalOneMinusLambda => "1 / (1 - lambda)",
        LegendreOrbitElementKind::LambdaMinusOneOverLambda => "(lambda - 1) / lambda",
        LegendreOrbitElementKind::LambdaOverLambdaMinusOne => "lambda / (lambda - 1)",
    }
}

pub(crate) fn format_legendre_parameter_conditioning(
    conditioning: LegendreParameterConditioning,
) -> &'static str {
    match conditioning {
        LegendreParameterConditioning::Generic => "generic",
        LegendreParameterConditioning::NearZero => "near zero",
        LegendreParameterConditioning::NearOne => "near one",
        LegendreParameterConditioning::NearInfinity => "near infinity",
    }
}

pub(crate) fn describe_invariant_recovery_interpretation(
    interpretation: InvariantRecoveryInterpretation,
) -> &'static str {
    match interpretation {
        InvariantRecoveryInterpretation::DirectAgreement => "direct agreement of g₂, g₃, Δ, and j",
        InvariantRecoveryInterpretation::SameModularClassButScaleSensitiveMismatch => {
            "same modular class via j, but scale-sensitive mismatch in g₂, g₃, or Δ"
        }
        InvariantRecoveryInterpretation::Inconsistent => {
            "inconsistent even at the modular-invariant level j"
        }
    }
}

pub(crate) fn format_fundamental_domain_step_reason(
    reason: FundamentalDomainReductionStepReason,
) -> &'static str {
    match reason {
        FundamentalDomainReductionStepReason::RealPartOutsideCenteredStrip => {
            "real part lay outside the centered strip"
        }
        FundamentalDomainReductionStepReason::NormLessThanOne => "point lay inside the unit circle",
    }
}

pub(crate) fn format_fundamental_domain_status(
    status: FundamentalDomainReductionStatus,
) -> &'static str {
    match status {
        FundamentalDomainReductionStatus::AlreadyReduced => "already reduced",
        FundamentalDomainReductionStatus::Reduced => "reduced",
        FundamentalDomainReductionStatus::StepLimitReached => "step limit reached",
    }
}

pub(crate) fn format_division_polynomial_status(
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

pub(crate) fn format_even_branch(branch: &EvenDivisionPolynomialVanishingBranch) -> &'static str {
    match branch {
        EvenDivisionPolynomialVanishingBranch::YApproxZero => "y(P) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::XCriterionApproxZero => "ε_n(x(P)) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::BothBranches => "both y(P) ≈ 0 and ε_n(x(P)) ≈ 0",
        EvenDivisionPolynomialVanishingBranch::NeitherBranch => {
            "neither y(P) nor ε_n(x(P)) is approximately zero"
        }
    }
}

pub(crate) fn format_analytic_cubic_model(curve: &AnalyticWeierstrassCurve) -> String {
    let mut equation = "y^2 = 4x^3".to_string();
    append_polynomial_term(&mut equation, -*curve.g2(), "x");
    append_polynomial_term(&mut equation, -*curve.g3(), "");
    equation
}
