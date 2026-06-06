use crate::elliptic_curves::analytic::{
    AnalyticDivisionPolynomialComparisonCase, AnalyticEvenDivisionPolynomialReport,
    AnalyticOddDivisionPolynomialReport, AnalyticTorsionPointApprox,
};
use crate::visualization::traits::Visualizable;

use crate::visualization::elliptic_curves::analytic::formatting::{
    format_complex_scalar_compact, format_division_polynomial_status, format_even_branch,
};
use crate::visualization::elliptic_curves::short_weierstrass::format_point_compact;

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
        format!("z = {}", format_complex_scalar_compact(point.torus_point().z())),
        format!("curve point = {}", format_point_compact(point.curve_point())),
        format!(
            "lies on curve under tolerance = {}",
            if point.lies_on_curve() { "yes" } else { "no" }
        ),
    ]
    .join("\n")
}

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
        format!("x = ℘(z) ≈ {}", format_complex_scalar_compact(report.x_value())),
        format!(
            "ψ_n(x) ≈ {}",
            format_complex_scalar_compact(report.psi_n_x())
        ),
        format!("|ψ_n(x)| = {:.6e}", report.absolute_value()),
        format!("status = {}", format_division_polynomial_status(report.status())),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

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
        format!("x = ℘(z) ≈ {}", format_complex_scalar_compact(report.x_value())),
        format!(
            "ε_n(x) ≈ {}",
            format_complex_scalar_compact(report.epsilon_n_x())
        ),
        format!("|ε_n(x)| = {:.6e}", report.absolute_value()),
        format!("branch = {}", format_even_branch(report.branch())),
        format!("status = {}", format_division_polynomial_status(report.status())),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            report.tolerance().absolute,
            report.tolerance().relative
        ),
    ]
    .join("\n")
}

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
        format!("ε_n(x) ≈ {}", format_complex_scalar_compact(self.epsilon_n_x()))
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
