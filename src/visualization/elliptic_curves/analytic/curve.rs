use crate::ComplexApprox;
use crate::elliptic_curves::analytic::{
    AnalyticCurveMembershipReport, AnalyticInvariants, AnalyticWeierstrassCurve,
    EllipticFunctionApproximation, HasPoleDistance, TorusToCurveMapResult, TorusToCurveValues,
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
    WeierstrassPApprox, WeierstrassPDerivativeApprox,
};
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::visualization::traits::Visualizable;

use crate::visualization::elliptic_curves::analytic::formatting::{
    append_polynomial_term, format_analytic_cubic_model as format_analytic_cubic_model_shared,
    format_complex_scalar_compact,
};
use crate::visualization::elliptic_curves::short_weierstrass::format_point_compact;

/// Formats the analytic cubic model `y² = 4x³ - g₂x - g₃` while suppressing
/// numerically negligible coefficients.
pub fn format_analytic_cubic_model(curve: &AnalyticWeierstrassCurve) -> String {
    format_analytic_cubic_model_shared(curve)
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

pub fn describe_weierstrass_p_approx(approximation: &WeierstrassPApprox) -> String {
    describe_elliptic_function_approximation(
        "Weierstrass ℘ approximation",
        approximation,
        Some(approximation.pole_distance()),
        false,
    )
}

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

pub fn describe_weierstrass_differential_equation(
    report: &WeierstrassDifferentialEquationReport,
) -> String {
    let mut lines = vec![
        "Weierstrass differential equation".to_string(),
        format!("z = {}", format_complex_scalar_compact(report.z())),
    ];

    match report.values() {
        TorusToCurveValues::Pole => lines.push("values = Pole".to_string()),
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
