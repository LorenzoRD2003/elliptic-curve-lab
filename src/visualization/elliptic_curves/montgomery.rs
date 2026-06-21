use core::fmt;

use crate::elliptic_curves::{MontgomeryCurve, traits::CurveModelConversion};
use crate::fields::traits::Field;
use crate::visualization::{
    elliptic_curves::{
        general_weierstrass::format_general_weierstrass_curve,
        short_weierstrass::format_curve as format_short_curve,
    },
    fields::traits::VisualizableField,
    traits::Visualizable,
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
    use crate::fields::{Fp, traits::Field};
    use crate::visualization::elliptic_curves::montgomery::{
        describe_montgomery_curve, describe_montgomery_general_embedding,
        describe_montgomery_short_reduction, format_montgomery_curve,
    };
    use crate::visualization::traits::Visualizable;

    type F3 = Fp<3>;
    type F5 = Fp<5>;

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
}
