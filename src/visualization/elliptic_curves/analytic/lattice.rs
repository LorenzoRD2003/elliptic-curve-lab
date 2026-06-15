use crate::elliptic_curves::analytic::{
    AnalyticInvariants, ComplexLattice, EisensteinSumApprox, ModularQParameter,
    TruncationConvergenceReport,
};
use crate::visualization::fields::complex_approx::format_complex;
use crate::visualization::traits::Visualizable;

use crate::visualization::elliptic_curves::analytic::formatting::format_complex_scalar_compact;

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
        format!("weight k = {}", sum.weight()),
        format!("truncation radius = {}", sum.truncation().radius()),
        format!("terms used = {}", sum.terms_used()),
        format!("value ≈ {}", format_complex_scalar_compact(sum.value())),
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

/// Describes one convergence-style comparison between two lattice truncations.
pub fn describe_truncation_convergence(report: &TruncationConvergenceReport) -> String {
    [
        "Truncation convergence".to_string(),
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
            self.weight(),
            format_complex_scalar_compact(self.value())
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
            format_complex_scalar_compact(self.g2()),
            format_complex_scalar_compact(self.g3()),
            format_complex_scalar_compact(self.j_invariant())
        )
    }

    fn describe(&self) -> String {
        crate::visualization::elliptic_curves::describe_analytic_invariants(self)
    }
}
