use num_complex::Complex64;

use crate::ComplexApprox;
use crate::elliptic_curves::{
    AnalyticCurveMembershipReport, AnalyticInvariants, AnalyticWeierstrassCurve, ComplexLattice,
    EisensteinSumApprox, EllipticFunctionApproximation, HasPoleDistance, ShortWeierstrassCurve,
    TorusToCurveMapResult, TorusToCurveValues, TruncationConvergenceReport,
    WeierstrassDifferentialEquationReport, WeierstrassDifferentialEquationStatus,
    WeierstrassPApprox, WeierstrassPDerivativeApprox,
};
use crate::visualization::Visualizable;
use crate::visualization::elliptic_curves::format_point_compact;
use crate::visualization::fields::format_complex;

fn is_small_real(value: f64) -> bool {
    value.abs() <= 1.0e-12
}

fn is_small_complex(value: &Complex64) -> bool {
    value.norm() <= 1.0e-12
}

fn format_complex_scalar_compact(value: &Complex64) -> String {
    if is_small_complex(value) {
        return "0".to_string();
    }

    if is_small_real(value.im) {
        return format!("{:.6}", value.re);
    }

    if is_small_real(value.re) {
        return format!("{:.6}i", value.im);
    }

    format_complex(value)
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

/// Describes a side-by-side comparison between two Eisenstein truncations.
pub fn describe_truncation_convergence(report: &TruncationConvergenceReport) -> String {
    [
        "Truncation comparison".to_string(),
        format!("small radius = {}", report.small.radius()),
        format!("large radius = {}", report.large.radius()),
        format!(
            "small value ≈ {}",
            format_complex_scalar_compact(&report.small_value)
        ),
        format!(
            "large value ≈ {}",
            format_complex_scalar_compact(&report.large_value)
        ),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(&report.difference)
        ),
        format!("|difference| = {:.6e}", report.absolute_difference),
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
        format!("point: {}", format_point_compact(&report.point)),
        format!("lhs ≈ {}", format_complex_scalar_compact(&report.lhs)),
        format!("rhs ≈ {}", format_complex_scalar_compact(&report.rhs)),
        format!(
            "difference ≈ {}",
            format_complex_scalar_compact(&report.difference)
        ),
        format!("|difference| = {:.6e}", report.absolute_error),
        format!(
            "holds under tolerance = {}",
            if report.is_on_curve { "yes" } else { "no" }
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

impl Visualizable for TruncationConvergenceReport {
    fn format_compact(&self) -> String {
        format!(
            "Δ_trunc ≈ {}",
            format_complex_scalar_compact(&self.difference)
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

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{
        describe_analytic_invariants, describe_complex_lattice, describe_eisenstein_sum,
        describe_torus_to_curve_map, describe_weierstrass_differential_equation,
        describe_weierstrass_p_approx, format_analytic_cubic_model,
        format_short_weierstrass_over_complex,
    };
    use crate::elliptic_curves::{
        AnalyticCurvePoint, AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice,
        EllipticFunctionTruncation, LatticeSumTruncation, UpperHalfPlanePoint, analytic_invariants,
        g4_sum, map_torus_point_to_curve, verify_weierstrass_differential_equation, weierstrass_p,
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
        assert!(map.describe().contains("Torus to curve map"));
        assert!(
            report
                .describe()
                .contains("Weierstrass differential equation")
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
}
