use core::fmt;
use std::hash::Hash;

use crate::elliptic_curves::short_weierstrass::isogenies::frobenius::FrobeniusVerschiebungFactorizationReport;
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::traits::Isogeny;
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::short_weierstrass::format_curve,
    isogenies::function_field_maps::{
        describe_differential_pullback_report, format_short_weierstrass_function_field_map,
    },
};

/// Describes the factorization `[p] = V \circ Frob_p` on one curve.
fn describe_frobenius_verschiebung_factorization_report<F>(
    report: &FrobeniusVerschiebungFactorizationReport<F>,
) -> String
where
    F: FiniteField + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash + PartialEq,
{
    [
        "Frobenius-Verschiebung factorization".to_string(),
        format!("curve: {}", format_curve(report.curve())),
        format!("characteristic p: {}", F::characteristic()),
        format!("scalar map: [{}]", report.scalar()),
        format!(
            "absolute Frobenius: {} -> {}",
            format_curve(report.frobenius().domain()),
            format_curve(report.frobenius().codomain())
        ),
        format!(
            "Verschiebung: {} -> {}",
            format_curve(report.verschiebung().domain_curve()),
            format_curve(report.verschiebung().codomain_curve())
        ),
        format!(
            "[p]^*(x), [p]^*(y): {}",
            format_short_weierstrass_function_field_map(report.multiplication_by_p_pullback())
        ),
        "certificate status: both duality relations verified".to_string(),
    ]
    .join("\n")
}

/// Explains the full factorization report in a more mathematical style.
fn explain_frobenius_verschiebung_factorization_report<F>(
    report: &FrobeniusVerschiebungFactorizationReport<F>,
) -> String
where
    F: FiniteField + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash + PartialEq,
{
    [
        "Characteristic-p factorization".to_string(),
        format!("Start with E = {}.", format_curve(report.curve())),
        format!(
            "The direct scalar pullback [{}]^* is computed from the generic point of E.",
            report.scalar()
        ),
        format!(
            "Absolute Frobenius gives Frob_p : {} -> {}.",
            format_curve(report.frobenius().domain()),
            format_curve(report.frobenius().codomain())
        ),
        "Verschiebung is reconstructed by inverting the absolute-Frobenius pullback on the two direct pullback coordinates of [p]^*.".to_string(),
        format!(
            "The stored Verschiebung pullback is {}.",
            format_short_weierstrass_function_field_map(report.verschiebung().as_function_field_map())
        ),
        "The certificate validates both identities: V ∘ Frob_p = [p]_E and Frob_p ∘ V = [p]_{E^(p)}.".to_string(),
        "Differential of Frobenius:".to_string(),
        describe_differential_pullback_report(
            &report
                .frobenius_differential_report()
                .expect("stored Frobenius in the factorization report should admit a differential report"),
        ),
        "Differential of Verschiebung:".to_string(),
        describe_differential_pullback_report(
            &report
                .verschiebung_differential_report()
                .expect("stored Verschiebung in the factorization report should admit a differential report"),
        ),
    ]
    .join("\n")
}

impl<F> Visualizable for FrobeniusVerschiebungFactorizationReport<F>
where
    F: FiniteField + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash + PartialEq,
{
    fn format_compact(&self) -> String {
        format!(
            "[{}] = V o Frob_p on {}",
            self.scalar(),
            format_curve(self.curve())
        )
    }

    fn describe(&self) -> String {
        describe_frobenius_verschiebung_factorization_report(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::isogenies::scalar_multiplication::ScalarMultiplicationIsogeny;
    use crate::visualization::Visualizable;

    type F5 = crate::fields::Fp5;

    fn curve() -> ShortWeierstrassCurve<F5> {
        ShortWeierstrassCurve::new(F5::from_i64(1), F5::from_i64(1)).expect("valid curve")
    }

    #[test]
    fn description_mentions_curve_frobenius_and_verschiebung() {
        let scalar = ScalarMultiplicationIsogeny::new(curve(), 5).expect("scalar should build");
        let report = scalar
            .frobenius_verschiebung_factorization_report()
            .expect("report should build");
        let description = describe_frobenius_verschiebung_factorization_report(&report);

        assert!(description.contains("Frobenius-Verschiebung factorization"));
        assert!(description.contains("absolute Frobenius"));
        assert!(description.contains("Verschiebung"));
        assert!(description.contains("certificate status"));
    }

    #[test]
    fn explanation_mentions_pth_root_story_and_differentials() {
        let scalar = ScalarMultiplicationIsogeny::new(curve(), 5).expect("scalar should build");
        let report = scalar
            .frobenius_verschiebung_factorization_report()
            .expect("report should build");
        let explanation = explain_frobenius_verschiebung_factorization_report(&report);

        assert!(explanation.contains("inverting the absolute-Frobenius pullback"));
        assert!(explanation.contains("V ∘ Frob_p = [p]_E"));
        assert!(explanation.contains("Differential of Frobenius"));
        assert!(explanation.contains("Differential of Verschiebung"));
    }

    #[test]
    fn visualizable_trait_reuses_compact_summary() {
        let scalar = ScalarMultiplicationIsogeny::new(curve(), 5).expect("scalar should build");
        let report = scalar
            .frobenius_verschiebung_factorization_report()
            .expect("report should build");

        assert!(report.format_compact().contains("[5] = V o Frob_p"));
        assert!(
            report
                .describe()
                .contains("Frobenius-Verschiebung factorization")
        );
    }
}
