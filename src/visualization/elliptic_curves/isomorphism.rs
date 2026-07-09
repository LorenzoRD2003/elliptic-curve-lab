use crate::visualization::*;
use core::fmt;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::{
        CurveIsomorphismError, ShortWeierstrassIsomorphism, ShortWeierstrassQuadraticTwist,
        TwistKind,
    },
    traits::CurveIsomorphism,
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::visualization::{
    VisualizableField,
    elliptic_curves::short_weierstrass::format_curve,
    shared::{format_field_elem as format_elem, yes_no},
    traits::Visualizable,
};

fn scaling_powers<F: Field>(u: &F::Elem) -> (F::Elem, F::Elem, F::Elem, F::Elem) {
    let u2 = F::square(u);
    let u3 = F::mul(&u2, u);
    let u4 = F::square(&u2);
    let u6 = F::mul(&u4, &u2);
    (u2, u3, u4, u6)
}

fn copy_curve<F: Field>(curve: &ShortWeierstrassCurve<F>) -> ShortWeierstrassCurve<F> {
    ShortWeierstrassCurve::new(curve.a().clone(), curve.b().clone())
        .expect("validated short-Weierstrass curves should stay valid when copied")
}

/// Formats a short-Weierstrass scaling isomorphism compactly.
fn format_isomorphism<F: Field>(isomorphism: &ShortWeierstrassIsomorphism<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    format!(
        "phi_u: {} -> {} with u = {}",
        format_curve(isomorphism.domain()),
        format_curve(isomorphism.codomain()),
        format_elem::<F>(isomorphism.scaling_factor())
    )
}

impl<F: Field> Visualizable for ShortWeierstrassIsomorphism<F>
where
    F::Elem: VisualizableField + fmt::Display,
{
    fn format_compact(&self) -> String {
        format_isomorphism(self)
    }

    fn describe(&self) -> String {
        describe_isomorphism(self)
    }
}

/// Describes a short-Weierstrass scaling isomorphism and its coefficient transport.
fn describe_isomorphism<F: Field>(isomorphism: &ShortWeierstrassIsomorphism<F>) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    let (u2, u3, u4, u6) = scaling_powers::<F>(isomorphism.scaling_factor());
    let codomain = isomorphism.codomain();

    [
        "Short-Weierstrass isomorphism".to_string(),
        format!("domain: {}", format_curve(isomorphism.domain())),
        format!("codomain: {}", format_curve(codomain)),
        format!("u: {}", format_elem::<F>(isomorphism.scaling_factor())),
        "map on affine points: (x, y) -> (u^2 x, u^3 y)".to_string(),
        format!("u^2: {}", format_elem::<F>(&u2)),
        format!("u^3: {}", format_elem::<F>(&u3)),
        format!("u^4: {}", format_elem::<F>(&u4)),
        format!("u^6: {}", format_elem::<F>(&u6)),
        format!(
            "coefficient transport: a' = u^4 a = {}, b' = u^6 b = {}",
            format_elem::<F>(codomain.a()),
            format_elem::<F>(codomain.b())
        ),
    ]
    .join("\n")
}

/// Explains the coefficient and coordinate scaling determined by `u`.
fn explain_short_weierstrass_scaling<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    u: &F::Elem,
) -> Result<String, CurveIsomorphismError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let isomorphism = ShortWeierstrassIsomorphism::new(copy_curve(curve), u.clone())?;
    let codomain = isomorphism.codomain();
    let (u2, u3, u4, u6) = scaling_powers::<F>(u);

    Ok([
        "Short-Weierstrass scaling".to_string(),
        format!("domain curve: {}", format_curve(curve)),
        format!("scaling parameter u: {}", format_elem::<F>(u)),
        "coordinate change: (x, y) -> (u^2 x, u^3 y)".to_string(),
        format!("u^2 = {}", format_elem::<F>(&u2)),
        format!("u^3 = {}", format_elem::<F>(&u3)),
        format!("u^4 = {}", format_elem::<F>(&u4)),
        format!("u^6 = {}", format_elem::<F>(&u6)),
        format!(
            "new coefficients: a' = u^4 a = {}, b' = u^6 b = {}",
            format_elem::<F>(codomain.a()),
            format_elem::<F>(codomain.b())
        ),
        format!("codomain curve: {}", format_curve(codomain)),
        "interpretation: this is an isomorphism over the current base field because u is invertible".to_string(),
    ]
    .join("\n"))
}

/// Explains the quadratic twist determined by `d`.
fn explain_quadratic_twist<F: SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    d: &F::Elem,
) -> Result<String, CurveIsomorphismError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let package = ShortWeierstrassQuadraticTwist::new(copy_curve(curve), d.clone())?;
    let d2 = F::square(d);
    let d3 = F::mul(&d2, d);

    let mut lines = vec![
        "Quadratic twist".to_string(),
        format!("original curve: {}", format_curve(curve)),
        format!("twist factor d: {}", format_elem::<F>(d)),
        format!("d^2 = {}", format_elem::<F>(&d2)),
        format!("d^3 = {}", format_elem::<F>(&d3)),
        format!("twisted curve: {}", format_curve(package.twist())),
        format!(
            "j-invariant preserved: {}",
            yes_no(curve.has_same_j_invariant(package.twist()))
        ),
    ];

    match package.kind() {
        TwistKind::Trivial => {
            lines.push("twist kind over current base field: trivial".to_string());
            lines.push(
                "interpretation: d is a square in the base field, so the twist is already isomorphic to the original curve over that field"
                    .to_string(),
            );
        }
        TwistKind::Quadratic => {
            lines.push("twist kind over current base field: quadratic".to_string());
            lines.push(
                "interpretation: d is not a square in the base field, so the twist typically becomes isomorphic only after adjoining sqrt(d)"
                    .to_string(),
            );
        }
    }

    Ok(lines.join("\n"))
}

/// Summarizes the comparison between two short-Weierstrass curves over a small enumerable field.
fn summarize_curve_comparison<F: EnumerableFiniteField>(
    left: &ShortWeierstrassCurve<F>,
    right: &ShortWeierstrassCurve<F>,
) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    let same_j = left.has_same_j_invariant(right);
    let isomorphism = left.find_isomorphism_to(right);
    let isomorphic_over_base = isomorphism.is_some();

    let interpretation = if same_j && isomorphic_over_base {
        "E and E' are related by a short-Weierstrass scaling over the current base field."
    } else if same_j {
        "E and E' become isomorphic over an algebraic extension, but no base-field scaling witness was found."
    } else {
        "E and E' do not even become isomorphic over an algebraic closure because their j-invariants differ."
    };

    let mut lines = vec![
        "Curve comparison".to_string(),
        "================".to_string(),
        String::new(),
        format!("E:  {}", format_curve(left)),
        format!("E': {}", format_curve(right)),
        String::new(),
        format!("j(E)  = {}", format_elem::<F>(&left.j_invariant())),
        format!("j(E') = {}", format_elem::<F>(&right.j_invariant())),
        String::new(),
        format!("same j-invariant: {}", yes_no(same_j)),
        format!(
            "isomorphic over base field: {}",
            yes_no(isomorphic_over_base)
        ),
    ];

    if let Some(isomorphism) = isomorphism {
        lines.push(format!(
            "base-field witness: u = {}",
            format_elem::<F>(isomorphism.scaling_factor())
        ));
    }

    lines.push(format!("interpretation: {interpretation}"));
    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        short_weierstrass::isomorphisms::{CurveIsomorphismError, ShortWeierstrassIsomorphism},
    };
    use crate::fields::traits::Field;

    type F7 = crate::fields::Fp7;
    type F19 = crate::fields::Fp19;

    fn f7_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f19_curve() -> ShortWeierstrassCurve<F19> {
        ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn compact_isomorphism_format_mentions_domain_codomain_and_u() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let formatted = format_isomorphism(&isomorphism);

        assert!(formatted.contains("phi_u"));
        assert!(formatted.contains("with u = 3"));
        assert!(formatted.contains("y^2 = x^3"));
    }

    #[test]
    fn isomorphism_description_mentions_map_and_coefficient_transport() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let description = describe_isomorphism(&isomorphism);

        assert!(description.contains("map on affine points"));
        assert!(description.contains("(x, y) -> (u^2 x, u^3 y)"));
        assert!(description.contains("coefficient transport"));
        assert!(description.contains("u^4"));
        assert!(description.contains("u^6"));
    }

    #[test]
    fn scaling_explanation_shows_powers_and_codomain() {
        let explanation = explain_short_weierstrass_scaling(&f7_curve(), &F7::from_i64(3))
            .expect("valid scaling should explain");

        assert!(explanation.contains("Short-Weierstrass scaling"));
        assert!(explanation.contains("u^2 = 2"));
        assert!(explanation.contains("u^3 = 6"));
        assert!(explanation.contains("codomain curve"));
    }

    #[test]
    fn scaling_explanation_rejects_noninvertible_u() {
        assert!(matches!(
            explain_short_weierstrass_scaling(&f7_curve(), &F7::zero()),
            Err(CurveIsomorphismError::NonInvertibleScale)
        ));
    }

    #[test]
    fn quadratic_twist_explanation_distinguishes_trivial_twists() {
        let explanation = explain_quadratic_twist(&f19_curve(), &F19::from_i64(4))
            .expect("square twist factor should explain");

        assert!(explanation.contains("Quadratic twist"));
        assert!(explanation.contains("j-invariant preserved: yes"));
        assert!(explanation.contains("twist kind over current base field: trivial"));
    }

    #[test]
    fn quadratic_twist_explanation_distinguishes_genuinely_quadratic_twists() {
        let explanation = explain_quadratic_twist(&f19_curve(), &F19::from_i64(2))
            .expect("non-square twist factor should explain");

        assert!(explanation.contains("twist kind over current base field: quadratic"));
        assert!(explanation.contains("adjoining sqrt(d)"));
    }

    #[test]
    fn curve_comparison_reports_same_j_but_not_base_field_isomorphic() {
        let left = f19_curve();
        let right = left
            .quadratic_twist(F19::from_i64(2))
            .expect("non-square twist should exist");
        let summary = summarize_curve_comparison(&left, &right);

        assert!(summary.contains("same j-invariant: yes"));
        assert!(summary.contains("isomorphic over base field: no"));
        assert!(summary.contains("become isomorphic over an algebraic extension"));
    }

    #[test]
    fn curve_comparison_reports_base_field_witness_when_present() {
        let left = f7_curve();
        let right = left
            .scaled_by(F7::from_i64(3))
            .expect("invertible scaling should produce a valid curve");
        let summary = summarize_curve_comparison(&left, &right);

        assert!(summary.contains("same j-invariant: yes"));
        assert!(summary.contains("isomorphic over base field: yes"));
        assert!(summary.contains("base-field witness: u = 3"));
    }

    #[test]
    fn curve_comparison_reports_different_j_when_curves_differ_over_closure() {
        let left = f7_curve();
        let right = ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::from_i64(1))
            .expect("valid curve");
        let summary = summarize_curve_comparison(&left, &right);

        assert!(summary.contains("same j-invariant: no"));
        assert!(summary.contains("do not even become isomorphic over an algebraic closure"));
    }
}
