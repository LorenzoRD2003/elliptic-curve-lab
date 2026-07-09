use core::fmt;

use crate::elliptic_curves::{GeneralWeierstrassCurve, traits::CurveModelConversion};
use crate::fields::traits::Field;
use crate::visualization::{
    Visualizable, VisualizableField,
    elliptic_curves::short_weierstrass::format_curve as format_short_curve,
    shared::{format_field_elem as format_elem, parenthesize_if_needed, yes_no},
};

fn compact_linear_term<F: Field>(coefficient: &F::Elem, variable: &str) -> Option<String>
where
    F::Elem: VisualizableField,
{
    if F::is_zero(coefficient) {
        None
    } else if F::eq(coefficient, &F::one()) {
        Some(variable.to_string())
    } else {
        Some(format!(
            "{}{}",
            parenthesize_if_needed(&format_elem::<F>(coefficient)),
            variable
        ))
    }
}

fn compact_constant_term<F: Field>(coefficient: &F::Elem) -> Option<String>
where
    F::Elem: VisualizableField,
{
    (!F::is_zero(coefficient)).then(|| parenthesize_if_needed(&format_elem::<F>(coefficient)))
}

/// Formats a general Weierstrass curve compactly.
pub(crate) fn format_general_weierstrass_curve<F: Field>(
    curve: &GeneralWeierstrassCurve<F>,
) -> String
where
    F::Elem: VisualizableField,
{
    let mut left_terms = vec!["y^2".to_string()];
    let mut right_terms = vec!["x^3".to_string()];

    if let Some(term) = compact_linear_term::<F>(curve.a1(), "xy") {
        left_terms.push(term);
    }
    if let Some(term) = compact_linear_term::<F>(curve.a3(), "y") {
        left_terms.push(term);
    }
    if let Some(term) = compact_linear_term::<F>(curve.a2(), "x^2") {
        right_terms.push(term);
    }
    if let Some(term) = compact_linear_term::<F>(curve.a4(), "x") {
        right_terms.push(term);
    }
    if let Some(term) = compact_constant_term::<F>(curve.a6()) {
        right_terms.push(term);
    }

    format!("{} = {}", left_terms.join(" + "), right_terms.join(" + "))
}

/// Describes a general Weierstrass curve in its native `a1,a2,a3,a4,a6`
/// presentation together with the classical invariants derived from it.
fn describe_general_weierstrass_curve<F: Field>(curve: &GeneralWeierstrassCurve<F>) -> String
where
    F::Elem: VisualizableField,
{
    [
        "General-Weierstrass curve".to_string(),
        format!("equation: {}", format_general_weierstrass_curve(curve)),
        format!("characteristic: {}", F::characteristic()),
        format!("a1: {}", format_elem::<F>(curve.a1())),
        format!("a2: {}", format_elem::<F>(curve.a2())),
        format!("a3: {}", format_elem::<F>(curve.a3())),
        format!("a4: {}", format_elem::<F>(curve.a4())),
        format!("a6: {}", format_elem::<F>(curve.a6())),
        format!("discriminant: {}", format_elem::<F>(&curve.discriminant())),
        format!("b2: {}", format_elem::<F>(&curve.b2())),
        format!("b4: {}", format_elem::<F>(&curve.b4())),
        format!("b6: {}", format_elem::<F>(&curve.b6())),
        format!("b8: {}", format_elem::<F>(&curve.b8())),
        format!("c4: {}", format_elem::<F>(&curve.c4())),
        format!("c6: {}", format_elem::<F>(&curve.c6())),
        format!("j-invariant: {}", format_elem::<F>(&curve.j_invariant())),
    ]
    .join("\n")
}

/// Describes the current explicit reduction route from the general model to a
/// short-Weierstrass companion.
fn describe_general_weierstrass_short_reduction<F: Field>(
    curve: &GeneralWeierstrassCurve<F>,
) -> String
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    let mut lines = vec![
        "General-to-short companion reduction".to_string(),
        format!("source curve: {}", format_general_weierstrass_curve(curve)),
    ];

    match curve.conversion_to_short_weierstrass() {
        Ok(conversion) => {
            lines.push("status: available in this characteristic".to_string());
            lines.push(format!(
                "target curve: {}",
                format_short_curve(conversion.target())
            ));
            lines.push(
                "route: algorithms that are still short-specific may delegate through this explicit companion"
                    .to_string(),
            );
            lines.push(format!(
                "invariants preserved: c4={}, c6={}, discriminant={}, j={}",
                yes_no(F::eq(&curve.c4(), &conversion.target().c4())),
                yes_no(F::eq(&curve.c6(), &conversion.target().c6())),
                yes_no(F::eq(
                    &curve.discriminant(),
                    &conversion.target().discriminant()
                )),
                yes_no(F::eq(
                    &curve.j_invariant(),
                    &conversion.target().j_invariant()
                )),
            ));
            lines.push("point transport: explicit in both directions".to_string());
        }
        Err(error) => {
            lines.push("status: unavailable in this characteristic".to_string());
            lines.push(format!("reason: {error}"));
            lines.push(
                "note: the general model remains valid here; only the classical short companion is unavailable"
                    .to_string(),
            );
        }
    }

    lines.join("\n")
}

impl<F: Field> Visualizable for GeneralWeierstrassCurve<F>
where
    F::Elem: VisualizableField + fmt::Display + Clone,
{
    fn format_compact(&self) -> String {
        format_general_weierstrass_curve(self)
    }

    fn describe(&self) -> String {
        describe_general_weierstrass_curve(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::GeneralWeierstrassCurve;
    use crate::fields::traits::Field;
    use crate::visualization::Visualizable;

    type F2 = crate::fields::Fp2;
    type F5 = crate::fields::Fp5;

    #[test]
    fn compact_formatter_shows_the_general_equation_terms() {
        let curve = GeneralWeierstrassCurve::<F5>::new(
            F5::one(),
            F5::one(),
            F5::one(),
            F5::one(),
            F5::zero(),
        )
        .expect("non-singular curve");

        let formatted = format_general_weierstrass_curve(&curve);

        assert!(formatted.contains("y^2 + xy + y = x^3 + x^2 + x"));
    }

    #[test]
    fn rich_description_mentions_native_general_invariants() {
        let curve = GeneralWeierstrassCurve::<F5>::new(
            F5::one(),
            F5::one(),
            F5::one(),
            F5::one(),
            F5::zero(),
        )
        .expect("non-singular curve");
        let description = describe_general_weierstrass_curve(&curve);

        assert!(description.contains("General-Weierstrass curve"));
        assert!(description.contains("a1:"));
        assert!(description.contains("b2:"));
        assert!(description.contains("j-invariant:"));
        assert_eq!(
            curve.format_compact(),
            format_general_weierstrass_curve(&curve)
        );
    }

    #[test]
    fn supported_reduction_description_mentions_the_short_companion() {
        let curve = GeneralWeierstrassCurve::<F5>::new(
            F5::one(),
            F5::one(),
            F5::one(),
            F5::one(),
            F5::zero(),
        )
        .expect("non-singular curve");
        let description = describe_general_weierstrass_short_reduction(&curve);

        assert!(description.contains("status: available"));
        assert!(description.contains("target curve: y^2 = x^3"));
        assert!(description.contains("invariants preserved:"));
        assert!(description.contains("point transport: explicit"));
    }

    #[test]
    fn unsupported_reduction_description_is_honest_about_characteristic_two() {
        let curve = GeneralWeierstrassCurve::<F2>::new(
            F2::one(),
            F2::zero(),
            F2::one(),
            F2::zero(),
            F2::one(),
        )
        .expect("non-singular curve");
        let description = describe_general_weierstrass_short_reduction(&curve);

        assert!(description.contains("status: unavailable"));
        assert!(description.contains("characteristic 2"));
        assert!(description.contains("general model remains valid"));
    }
}
