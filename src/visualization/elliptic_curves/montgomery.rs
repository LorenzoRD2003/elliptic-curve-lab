use crate::visualization::*;
use core::fmt;

use crate::elliptic_curves::montgomery::{
    MontgomeryLadderReport, MontgomeryXzPoint, NormalizedMontgomeryCurve,
};
use crate::elliptic_curves::{MontgomeryCurve, traits::CurveModelConversion};
use crate::fields::traits::SqrtField;
use crate::visualization::{
    Visualizable,
    elliptic_curves::{
        general_weierstrass::format_general_weierstrass_curve,
        short_weierstrass::format_curve as format_short_curve,
    },
};

fn format_elem<F: Field>(value: &F::Elem) -> String
where
    F::Elem: VisualizableField,
{
    value.format_elem()
}

fn parenthesize_if_needed(text: &str) -> String {
    if text.contains(' ') || text.starts_with('-') || text.contains('/') {
        format!("({text})")
    } else {
        text.to_string()
    }
}

/// Formats a Montgomery curve compactly.
pub fn format_montgomery_curve<F: Field>(curve: &MontgomeryCurve<F>) -> String
where
    F::Elem: VisualizableField,
{
    format!(
        "{}y^2 = x^3 + {}x^2 + x",
        parenthesize_if_needed(&format_elem::<F>(curve.b())),
        parenthesize_if_needed(&format_elem::<F>(curve.a())),
    )
}

/// Formats one projective Montgomery `x`-line value compactly.
pub fn format_montgomery_xz_point<F: Field>(point: &MontgomeryXzPoint<F>) -> String
where
    F::Elem: VisualizableField,
{
    match point {
        MontgomeryXzPoint::Infinity => "O_x".to_string(),
        MontgomeryXzPoint::Finite { x, z } => format!(
            "({} : {})",
            parenthesize_if_needed(&format_elem::<F>(x)),
            parenthesize_if_needed(&format_elem::<F>(z)),
        ),
    }
}

/// Describes a Montgomery curve in its native `A,B` presentation together
/// with the classical invariants derived from it.
pub fn describe_montgomery_curve<F: Field>(curve: &MontgomeryCurve<F>) -> String
where
    F::Elem: VisualizableField,
{
    [
        "Montgomery curve".to_string(),
        format!("equation: {}", format_montgomery_curve(curve)),
        format!("characteristic: {}", F::characteristic()),
        format!("A: {}", format_elem::<F>(curve.a())),
        format!("B: {}", format_elem::<F>(curve.b())),
        format!("discriminant: {}", format_elem::<F>(&curve.discriminant())),
        format!("c4: {}", format_elem::<F>(&curve.c4())),
        format!("c6: {}", format_elem::<F>(&curve.c6())),
        format!("j-invariant: {}", format_elem::<F>(&curve.j_invariant())),
    ]
    .join("\n")
}

/// Describes the current explicit reduction route from the Montgomery model to
/// a short-Weierstrass companion.
pub fn describe_montgomery_short_reduction<F: Field>(curve: &MontgomeryCurve<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let mut lines = vec![
        "Montgomery-to-short companion reduction".to_string(),
        format!("source curve: {}", format_montgomery_curve(curve)),
    ];

    match curve.conversion_to_short_weierstrass() {
        Ok(conversion) => {
            lines.push("status: available in this characteristic".to_string());
            lines.push(format!(
                "target curve: {}",
                format_short_curve(conversion.target())
            ));
            lines.push(
                "route: deeper finite-field algorithms may delegate through this explicit short companion"
                    .to_string(),
            );
            lines.push(format!(
                "invariants preserved: c4={}, c6={}, discriminant={}, j={}",
                if F::eq(&curve.c4(), &conversion.target().c4()) {
                    "yes"
                } else {
                    "no"
                },
                if F::eq(&curve.c6(), &conversion.target().c6()) {
                    "yes"
                } else {
                    "no"
                },
                if F::eq(&curve.discriminant(), &conversion.target().discriminant()) {
                    "yes"
                } else {
                    "no"
                },
                if F::eq(&curve.j_invariant(), &conversion.target().j_invariant()) {
                    "yes"
                } else {
                    "no"
                },
            ));
            lines.push("point transport: explicit in both directions".to_string());
        }
        Err(error) => {
            lines.push("status: unavailable in this characteristic".to_string());
            lines.push(format!("reason: {error}"));
            lines.push(
                "note: the Montgomery model itself remains valid here; only the classical short companion is unavailable"
                    .to_string(),
            );
        }
    }

    lines.join("\n")
}

/// Describes the direct inclusion of the Montgomery model into the general
/// Weierstrass family.
pub fn describe_montgomery_general_embedding<F: Field>(curve: &MontgomeryCurve<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let general = curve.as_general_weierstrass();

    [
        "Montgomery-to-general embedding".to_string(),
        format!("source curve: {}", format_montgomery_curve(curve)),
        format!(
            "target curve: {}",
            format_general_weierstrass_curve(&general)
        ),
        "route: direct affine rescaling, without passing through the short companion".to_string(),
    ]
    .join("\n")
}

/// Describes one Montgomery ladder report on the normalized model
///
/// `v^2 = x^3 + A x^2 + x`.
pub fn describe_normalized_montgomery_ladder_report<F: Field>(
    curve: &NormalizedMontgomeryCurve<F>,
    report: &MontgomeryLadderReport<F>,
) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    [
        "Normalized Montgomery ladder".to_string(),
        format!("curve: {}", curve),
        "route: native x-only Montgomery ladder on the normalized B = 1 model".to_string(),
        format!("input x(P): {}", format_elem::<F>(report.base_x())),
        format!("scalar n: {}", report.scalar()),
        format!("x([n]P): {}", format_montgomery_xz_point(report.multiple_x())),
        format!(
            "x([n+1]P): {}",
            format_montgomery_xz_point(report.next_multiple_x())
        ),
        "state invariant: the tracked pair is (x([n]P), x([n+1]P)) with fixed difference P at the point level".to_string(),
        "scope: this report is x-only; it does not determine the sign of y([n]P)".to_string(),
    ]
    .join("\n")
}

/// Describes one Montgomery ladder report from the source model
///
/// `B y^2 = x^3 + A x^2 + x`.
pub fn describe_montgomery_ladder_report<F: Field>(
    curve: &MontgomeryCurve<F>,
    report: &MontgomeryLadderReport<F>,
) -> String
where
    F: SqrtField,
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let mut lines = vec![
        "Montgomery ladder".to_string(),
        format!("source curve: {}", format_montgomery_curve(curve)),
        format!("input x(P): {}", format_elem::<F>(report.base_x())),
        format!("scalar n: {}", report.scalar()),
        format!(
            "x([n]P): {}",
            format_montgomery_xz_point(report.multiple_x())
        ),
        format!(
            "x([n+1]P): {}",
            format_montgomery_xz_point(report.next_multiple_x())
        ),
        "scope: the output is an x-coordinate class, not a signed affine point".to_string(),
    ];

    match curve.try_as_normalized_montgomery() {
        Ok(normalized) => {
            lines.push(
                "route: delegated through the same-field normalization to the native B = 1 Montgomery ladder"
                    .to_string(),
            );
            lines.push(format!("normalized companion: {}", normalized));
        }
        Err(error) => {
            lines.push(format!(
                "status: unavailable on this source curve ({error})"
            ));
        }
    }

    lines.join("\n")
}

impl<F> Visualizable for MontgomeryCurve<F>
where
    F: Field,
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    fn format_compact(&self) -> String {
        format_montgomery_curve(self)
    }

    fn describe(&self) -> String {
        describe_montgomery_curve(self)
    }
}

#[cfg(test)]
mod tests {

    use crate::elliptic_curves::MontgomeryCurve;
    use crate::fields::traits::Field;
    use crate::visualization::elliptic_curves::montgomery::{
        describe_montgomery_curve, describe_montgomery_general_embedding,
        describe_montgomery_ladder_report, describe_montgomery_short_reduction,
        describe_normalized_montgomery_ladder_report, format_montgomery_curve,
        format_montgomery_xz_point,
    };
    use crate::visualization::traits::Visualizable;

    type F3 = crate::fields::Fp3;
    type F5 = crate::fields::Fp5;

    fn f3_curve() -> MontgomeryCurve<F3> {
        MontgomeryCurve::<F3>::new(F3::zero(), F3::one()).expect("non-singular curve")
    }

    fn f5_curve() -> MontgomeryCurve<F5> {
        MontgomeryCurve::<F5>::new(F5::one(), F5::one()).expect("non-singular curve")
    }

    #[test]
    fn compact_formatter_shows_the_montgomery_equation_terms() {
        let curve = f5_curve();

        let formatted = format_montgomery_curve(&curve);

        assert!(formatted.contains("y^2 = x^3 + "));
        assert!(formatted.contains("x^2 + x"));
    }

    #[test]
    fn rich_description_mentions_native_montgomery_invariants() {
        let curve = f5_curve();
        let description = describe_montgomery_curve(&curve);

        assert!(description.contains("Montgomery curve"));
        assert!(description.contains("A:"));
        assert!(description.contains("B:"));
        assert!(description.contains("j-invariant:"));
        assert_eq!(curve.format_compact(), format_montgomery_curve(&curve));
    }

    #[test]
    fn supported_reduction_description_mentions_the_short_companion() {
        let curve = f5_curve();
        let description = describe_montgomery_short_reduction(&curve);

        assert!(description.contains("status: available"));
        assert!(description.contains("target curve: y^2 = x^3"));
        assert!(description.contains("invariants preserved:"));
        assert!(description.contains("point transport: explicit"));
    }

    #[test]
    fn unsupported_reduction_description_mentions_the_characteristic_limit() {
        let curve = f3_curve();
        let description = describe_montgomery_short_reduction(&curve);

        assert!(description.contains("status: unavailable"));
        assert!(description.contains("reason:"));
        assert!(description.contains("Montgomery model itself remains valid"));
    }

    #[test]
    fn general_embedding_description_mentions_the_direct_general_view() {
        let curve = f5_curve();
        let description = describe_montgomery_general_embedding(&curve);

        assert!(description.contains("Montgomery-to-general embedding"));
        assert!(description.contains("target curve: y^2 = x^3 + x^2 + x"));
        assert!(description.contains("direct affine rescaling"));
    }

    #[test]
    fn xz_point_formatter_uses_the_projective_x_line_story() {
        let curve = f5_curve();
        let report = curve
            .try_ladder_x_report(F5::from_i64(2), 3)
            .expect("ladder should be available in this example");

        let formatted = format_montgomery_xz_point(report.multiple_x());

        assert!(formatted.starts_with("(") || formatted == "O_x");
    }

    #[test]
    fn ladder_descriptions_surface_x_only_scope_and_route_honestly() {
        let curve = f5_curve();
        let report = curve
            .try_ladder_x_report(F5::from_i64(2), 3)
            .expect("ladder should be available in this example");
        let normalized = curve
            .try_as_normalized_montgomery()
            .expect("normalization should be available in this example");

        let source_description = describe_montgomery_ladder_report(&curve, &report);
        let normalized_description =
            describe_normalized_montgomery_ladder_report(&normalized, &report);

        assert!(source_description.contains("Montgomery ladder"));
        assert!(source_description.contains("input x(P):"));
        assert!(source_description.contains("x([n]P):"));
        assert!(source_description.contains("x-coordinate class"));
        assert!(source_description.contains("delegated through"));
        assert!(normalized_description.contains("Normalized Montgomery ladder"));
        assert!(normalized_description.contains("native x-only"));
        assert!(normalized_description.contains("state invariant:"));
        assert!(normalized_description.contains("x-only"));
    }
}
