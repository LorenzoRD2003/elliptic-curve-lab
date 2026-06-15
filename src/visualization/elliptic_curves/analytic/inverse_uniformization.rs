use crate::elliptic_curves::analytic::{
    InvariantRecoveryInterpretation, InvariantRecoveryValidationReport,
    InverseUniformizationJValidationReport, PointRoundTripValidationConfig,
    PointRoundTripValidationReport,
};
use crate::visualization::traits::Visualizable;

use crate::visualization::elliptic_curves::analytic::formatting::format_analytic_cubic_model;
use crate::visualization::elliptic_curves::analytic::formatting::{
    describe_invariant_recovery_interpretation, format_complex_scalar_compact,
    format_decimal_diagnostic,
};
use crate::visualization::elliptic_curves::short_weierstrass::format_point_compact;

pub fn describe_inverse_uniformization_j_validation_report(
    report: &InverseUniformizationJValidationReport,
) -> String {
    [
        "Inverse-uniformization j-validation report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("τ ≈ {}", format_complex_scalar_compact(report.tau().tau())),
        format!(
            "lattice basis ≈ ({}, {})",
            format_complex_scalar_compact(report.lattice().omega1()),
            format_complex_scalar_compact(report.lattice().omega2())
        ),
        format!(
            "recovered g₂ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().g2())
        ),
        format!(
            "recovered g₃ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().g3())
        ),
        format!(
            "recovered Δ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().discriminant())
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
            "lattice truncation radius = {}",
            report.lattice_truncation().radius()
        ),
        format!(
            "agrees under tolerance = {}",
            if report.agrees_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
        "this validates the modular j-class seen from the recovered τ, not the full scale-sensitive normalization".to_string(),
    ]
    .join("\n")
}

pub fn describe_invariant_recovery_validation_report(
    report: &InvariantRecoveryValidationReport,
) -> String {
    [
        "Invariant recovery validation report".to_string(),
        format!("curve = {}", format_analytic_cubic_model(report.curve())),
        format!("τ ≈ {}", format_complex_scalar_compact(report.tau().tau())),
        format!(
            "lattice basis ≈ ({}, {})",
            format_complex_scalar_compact(report.periods().omega1()),
            format_complex_scalar_compact(report.periods().omega2())
        ),
        format!(
            "recovered g₂ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().g2())
        ),
        format!(
            "curve-side g₂ ≈ {}",
            format_complex_scalar_compact(report.g2_comparison().right())
        ),
        format!(
            "Δg₂ ≈ {}",
            format_complex_scalar_compact(report.g2_comparison().difference())
        ),
        format!(
            "recovered g₃ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().g3())
        ),
        format!(
            "curve-side g₃ ≈ {}",
            format_complex_scalar_compact(report.g3_comparison().right())
        ),
        format!(
            "Δg₃ ≈ {}",
            format_complex_scalar_compact(report.g3_comparison().difference())
        ),
        format!(
            "recovered Δ ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().discriminant())
        ),
        format!(
            "curve-side Δ ≈ {}",
            format_complex_scalar_compact(report.discriminant_comparison().right())
        ),
        format!(
            "ΔΔ ≈ {}",
            format_complex_scalar_compact(report.discriminant_comparison().difference())
        ),
        format!(
            "recovered j ≈ {}",
            format_complex_scalar_compact(report.recovered_invariants().j_invariant())
        ),
        format!(
            "curve-side j ≈ {}",
            format_complex_scalar_compact(report.j_comparison().right())
        ),
        format!(
            "Δj ≈ {}",
            format_complex_scalar_compact(report.j_comparison().difference())
        ),
        format!(
            "interpretation = {}",
            describe_invariant_recovery_interpretation(report.interpretation())
        ),
        format!(
            "direct scale-sensitive agreement = {}",
            if report.direct_scale_sensitive_agreement() {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "same j-invariant approximately = {}",
            if report.same_j_invariant_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "lattice truncation radius = {}",
            report.lattice_truncation().radius()
        ),
        "g₂, g₃, and Δ are scale-sensitive, while j is homothety-invariant".to_string(),
    ]
    .join("\n")
}

pub fn describe_point_roundtrip_validation_config(
    config: &PointRoundTripValidationConfig,
) -> String {
    [
        "Point roundtrip validation config".to_string(),
        format!(
            "lattice truncation radius = {}",
            config.lattice_truncation().radius()
        ),
        format!(
            "elliptic-function truncation radius = {}",
            config.function_truncation().radius()
        ),
        format!(
            "tolerance = abs {:.3e}, rel {:.3e}",
            config.tolerance().absolute,
            config.tolerance().relative
        ),
        "this config tunes the forward check z -> (wp(z), wp'(z)), not the inverse Abel-Jacobi quadrature itself".to_string(),
    ]
    .join("\n")
}

pub fn describe_point_roundtrip_validation_report(
    report: &PointRoundTripValidationReport,
) -> String {
    let source_point = format_point_compact(report.point());
    let recovered_point = format_point_compact(report.recovered_curve_point());

    [
        "Point inverse-uniformization roundtrip report".to_string(),
        format!("source point P = {source_point}"),
        format!(
            "recovered torus representative z_P ≈ {}",
            format_complex_scalar_compact(report.reduced_representative())
        ),
        format!(
            "torus coordinates in the recovered basis ≈ ({}, {})",
            format_decimal_diagnostic(report.torus_point().coordinate().u()),
            format_decimal_diagnostic(report.torus_point().coordinate().v())
        ),
        format!("recovered curve point ≈ {recovered_point}"),
        format!("x residual norm = {:.6e}", report.x_residual_norm()),
        format!("y residual norm = {:.6e}", report.y_residual_norm()),
        format!(
            "lattice truncation radius = {}",
            report.lattice_truncation().radius()
        ),
        format!(
            "elliptic-function truncation radius = {}",
            report.function_truncation().radius()
        ),
        format!(
            "agrees under tolerance = {}",
            if report.agrees_approximately() {
                "yes"
            } else {
                "no"
            }
        ),
        "this report reuses the recovered torus class instead of rebuilding a second parallel inverse path".to_string(),
    ]
    .join("\n")
}

impl Visualizable for InverseUniformizationJValidationReport {
    fn format_compact(&self) -> String {
        format!(
            "Δj_inverse ≈ {}",
            format_complex_scalar_compact(self.difference())
        )
    }

    fn describe(&self) -> String {
        describe_inverse_uniformization_j_validation_report(self)
    }
}

impl Visualizable for InvariantRecoveryValidationReport {
    fn format_compact(&self) -> String {
        match self.interpretation() {
            InvariantRecoveryInterpretation::DirectAgreement => {
                "invariants: direct agreement".to_string()
            }
            InvariantRecoveryInterpretation::SameModularClassButScaleSensitiveMismatch => {
                format!(
                    "invariants: same j, Δg₂ ≈ {}, Δg₃ ≈ {}",
                    format_complex_scalar_compact(self.g2_comparison().difference()),
                    format_complex_scalar_compact(self.g3_comparison().difference())
                )
            }
            InvariantRecoveryInterpretation::Inconsistent => {
                format!(
                    "invariants: inconsistent, Δj ≈ {}",
                    format_complex_scalar_compact(self.j_comparison().difference())
                )
            }
        }
    }

    fn describe(&self) -> String {
        describe_invariant_recovery_validation_report(self)
    }
}

impl Visualizable for PointRoundTripValidationConfig {
    fn format_compact(&self) -> String {
        format!(
            "point-roundtrip config: r_Λ = {}, r_fun = {}",
            self.lattice_truncation().radius(),
            self.function_truncation().radius()
        )
    }

    fn describe(&self) -> String {
        describe_point_roundtrip_validation_config(self)
    }
}

impl Visualizable for PointRoundTripValidationReport {
    fn format_compact(&self) -> String {
        format!(
            "point roundtrip: {}",
            if self.agrees_approximately() {
                "agrees"
            } else {
                "does not agree"
            }
        )
    }

    fn describe(&self) -> String {
        describe_point_roundtrip_validation_report(self)
    }
}
