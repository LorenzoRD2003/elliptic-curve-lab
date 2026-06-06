use crate::elliptic_curves::analytic::{
    FundamentalDomainReductionReport, FundamentalDomainReductionStep, JInvariantComparisonReport,
    ModularInvarianceReport, ModularMatrix,
};
use crate::numerics::HasComplexApproxComparison;
use crate::visualization::Visualizable;

use crate::visualization::elliptic_curves::analytic::formatting::{
    format_complex_scalar_compact, format_fundamental_domain_status,
    format_fundamental_domain_step_reason,
};

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

/// Describes one modular matrix in `SL₂(ℤ)`.
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
        "determinant = 1".to_string(),
        "action on τ: (aτ + b) / (cτ + d)".to_string(),
    ]
    .join("\n")
}

/// Describes one modular-invariance comparison report.
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

/// Describes one reduction step toward the standard fundamental domain.
pub fn describe_fundamental_domain_reduction_step(step: &FundamentalDomainReductionStep) -> String {
    [
        "Fundamental-domain reduction step".to_string(),
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
        format!(
            "accumulated matrix = [[{}, {}], [{}, {}]]",
            step.applied_matrix().a(),
            step.applied_matrix().b(),
            step.applied_matrix().c(),
            step.applied_matrix().d()
        ),
    ]
    .join("\n")
}

/// Describes one full reduction report to the standard fundamental domain.
pub fn describe_fundamental_domain_reduction_report(
    report: &FundamentalDomainReductionReport,
) -> String {
    [
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
            "status = {}",
            format_fundamental_domain_status(report.status())
        ),
        format!("steps used = {}", report.steps().len()),
        format!(
            "accumulated matrix = [[{}, {}], [{}, {}]]",
            report.accumulated_matrix().a(),
            report.accumulated_matrix().b(),
            report.accumulated_matrix().c(),
            report.accumulated_matrix().d()
        ),
    ]
    .join("\n")
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
