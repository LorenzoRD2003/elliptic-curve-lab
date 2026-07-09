use core::fmt;
use std::hash::Hash;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isogenies::{DualVeluIsogeny, VeluIsogeny},
    traits::CurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, Field, FiniteField, SqrtField};
use crate::isogenies::{
    dual_report::{DualIsogenyReport, DualityKind},
    error::IsogenyError,
    scalar_multiplication::ScalarMultiplicationIsogeny,
    traits::Isogeny,
};
use crate::visualization::{
    Visualizable, VisualizableField, elliptic_curves::short_weierstrass::format_curve,
    shared::yes_no,
};

fn duality_kind_label(kind: DualityKind) -> &'static str {
    match kind {
        DualityKind::SeparableClassical => "separable classical",
        DualityKind::FrobeniusVerschiebung => "Frobenius/Verschiebung",
        DualityKind::MixedOrPartial => "mixed or partial",
    }
}

/// Describes an explicit composition between short-Weierstrass isogenies.
///
/// This helper is intentionally compact. It focuses on the three pieces that
/// matter most in the current implementation:
///
/// - the domain curve
/// - the codomain curve
/// - the multiplicative degree rule for the composed map
fn describe_composition<M, F>(composition: &M) -> String
where
    F: FiniteField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
    M: Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>,
{
    [
        "Composition of explicit isogenies".to_string(),
        format!("domain: {}", format_curve(composition.domain())),
        format!("codomain: {}", format_curve(composition.codomain())),
        format!("degree: {}", composition.degree()),
        "evaluation order: second(first(P))".to_string(),
    ]
    .join("\n")
}

/// Describes the scalar-multiplication self-isogeny `[n] : E -> E`.
fn describe_scalar_multiplication_isogeny<F>(
    map: &ScalarMultiplicationIsogeny<ShortWeierstrassCurve<F>>,
) -> String
where
    F: FiniteField + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    [
        "Scalar-multiplication isogeny".to_string(),
        format!("curve: {}", format_curve(map.domain())),
        format!("scalar: {}", map.scalar()),
        format!("degree: {}", map.degree()),
        format!("kernel size on E(F_q): {}", map.kernel_points().len()),
    ]
    .join("\n")
}

/// Describes the currently implemented short-Weierstrass dual-isogeny object.
///
/// The present dual search returns a map represented as:
///
/// - a Vélu isogeny built from a kernel on `E'(F_q)`
/// - a base-field isomorphism back to the original curve `E`
fn describe_dual_isogeny<F: Field + Clone>(dual: &DualVeluIsogeny<F>) -> String
where
    F::Elem: Clone + Eq + Hash,
{
    [
        format!("degree: {}", dual.degree()),
        "constructed by Velu from a kernel on E'".to_string(),
        "followed by an isomorphism back to E".to_string(),
    ]
    .join("\n")
}

impl<F: Field + Clone> Visualizable for DualVeluIsogeny<F>
where
    F::Elem: Clone + Eq + Hash,
{
    fn format_compact(&self) -> String {
        format!("dual Vélu isogeny of degree {}", self.degree())
    }

    fn describe(&self) -> String {
        describe_dual_isogeny(self)
    }
}

/// Describes a structured dual-isogeny report.
fn describe_dual_isogeny_report<Domain: CurveModel, Codomain: CurveModel>(
    report: &DualIsogenyReport<Domain, Codomain>,
) -> String {
    [
        "Dual isogeny report".to_string(),
        format!(
            "duality kind: {}",
            duality_kind_label(report.duality_kind())
        ),
        format!("deg(phi): {}", report.phi_degree()),
        format!("deg(phi_hat): {}", report.dual_degree()),
        format!("phi kernel: {}", report.phi_kernel_summary().short_label()),
        format!(
            "phi_hat kernel: {}",
            report.dual_kernel_summary().short_label()
        ),
    ]
    .join("\n")
}

impl<Domain: CurveModel, Codomain: CurveModel> Visualizable
    for DualIsogenyReport<Domain, Codomain>
{
    fn format_compact(&self) -> String {
        format!(
            "dual isogeny report: deg(phi) = {}, deg(phi_hat) = {}",
            self.phi_degree(),
            self.dual_degree()
        )
    }

    fn describe(&self) -> String {
        describe_dual_isogeny_report(self)
    }
}

/// Explains the expected duality relations in compact algebraic form.
fn explain_dual_relation<F: Field + Clone>(
    phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    dual: &DualVeluIsogeny<F>,
) -> String
where
    F::Elem: Clone + Eq + Hash,
{
    let degree = phi.degree();

    [
        "Dual relation".to_string(),
        "phi: E -> E'".to_string(),
        "phi_hat: E' -> E".to_string(),
        format!("degree(phi) = {}", degree),
        format!("degree(phi_hat) = {}", dual.degree()),
        format!("phi_hat o phi = [{}]_E", degree),
        format!("phi o phi_hat = [{}]_E'", degree),
    ]
    .join("\n")
}

/// Summarizes the exhaustive rational-point verification of a dual candidate.
///
/// This helper is intentionally small-scale and stage-specific: it checks
/// the left and right duality relations by full enumeration on the rational
/// points of `E(F_q)` and `E'(F_q)`.
fn summarize_dual_verification<F>(
    phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    dual: &DualVeluIsogeny<F>,
) -> Result<String, IsogenyError>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    let report = DualVeluIsogeny::dual_report(phi, dual)?;
    let degree = report.phi_degree();
    let left_holds = report.left_relation_holds();
    let right_holds = report.right_relation_holds();
    let composed_degree = report.phi_degree() * report.dual_degree();
    let scalar_degree = degree * degree;

    Ok([
        format!(
            "duality kind: {}",
            duality_kind_label(report.duality_kind())
        ),
        format!(
            "phi_hat(phi(Q)) = [{}]Q on all Q in E(Fp): {}",
            degree,
            yes_no(left_holds)
        ),
        format!(
            "phi(phi_hat(R)) = [{}]R on all R in E'(Fp): {}",
            degree,
            yes_no(right_holds)
        ),
        format!(
            "deg(phi_hat o phi) = {} = deg([{}]): {}",
            composed_degree,
            degree,
            yes_no(composed_degree == scalar_degree)
        ),
    ]
    .join("\n"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::{
        ShortWeierstrassCurve, short_weierstrass::isogenies::VeluIsogeny, traits::AffineCurveModel,
    };
    use crate::isogenies::{
        comparison::maps_equal_exhaustively, composition::ComposedIsogeny,
        scalar_multiplication::ScalarMultiplicationIsogeny, traits::Isogeny,
    };

    type F29 = crate::fields::Fp29;
    type F41 = crate::fields::Fp41;
    type Curve29 = ShortWeierstrassCurve<F29>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn curve_f29() -> Curve29 {
        Curve29::new(F29::from_i64(2), F29::from_i64(2)).expect("valid curve")
    }

    fn curve_f41() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn degree_three_phi() -> VeluIsogeny<Curve29> {
        let curve = curve_f29();
        let generator = curve
            .point(F29::from_i64(10), F29::from_i64(23))
            .expect("sample generator should lie on the curve");
        VeluIsogeny::from_generator(curve, generator).expect("sample Vélu isogeny should build")
    }

    fn strict_composition() -> ComposedIsogeny<
        VeluIsogeny<Curve41>,
        ScalarMultiplicationIsogeny<Curve41>,
        Curve41,
        Curve41,
        Curve41,
    > {
        let phi = {
            let curve = curve_f41();
            let generator = curve
                .point(F41::from_i64(40), F41::from_i64(0))
                .expect("sample generator should lie on the curve");
            VeluIsogeny::from_generator(curve, generator).expect("sample Vélu isogeny should build")
        };
        let identity_like = ScalarMultiplicationIsogeny::new(phi.codomain().clone(), 1)
            .expect("scalar-one isogeny should build");
        let composition = ComposedIsogeny::new_strict(phi.clone(), identity_like)
            .expect("composition should build");

        assert_eq!(
            maps_equal_exhaustively::<_, _, Curve41, Curve41>(&composition, &phi),
            Ok(true)
        );

        composition
    }

    #[test]
    fn composition_description_mentions_degree_and_evaluation_order() {
        let description = describe_composition(&strict_composition());

        assert!(description.contains("Composition of explicit isogenies"));
        assert!(description.contains("degree:"));
        assert!(description.contains("second(first(P))"));
    }

    #[test]
    fn scalar_multiplication_description_mentions_scalar_and_square_degree() {
        let description = describe_scalar_multiplication_isogeny(
            &ScalarMultiplicationIsogeny::new(curve_f41(), 3).expect("scalar map should build"),
        );

        assert!(description.contains("Scalar-multiplication isogeny"));
        assert!(description.contains("scalar: 3"));
        assert!(description.contains("degree: 9"));
    }

    #[test]
    fn dual_description_mentions_velu_then_isomorphism_construction() {
        let dual = degree_three_phi()
            .find_dual_exhaustively()
            .expect("dual should be found");
        let description = describe_dual_isogeny(&dual);

        assert!(description.contains("degree: 3"));
        assert!(description.contains("constructed by Velu from a kernel on E'"));
        assert!(description.contains("followed by an isomorphism back to E"));
    }

    #[test]
    fn dual_relation_explanation_mentions_both_duality_equations() {
        let phi = degree_three_phi();
        let dual = phi.find_dual_exhaustively().expect("dual should be found");
        let explanation = explain_dual_relation(&phi, &dual);

        assert!(explanation.contains("phi_hat o phi = [3]_E"));
        assert!(explanation.contains("phi o phi_hat = [3]_E'"));
    }

    #[test]
    fn dual_verification_summary_reports_yes_for_the_degree_three_example() {
        let phi = degree_three_phi();
        let dual = phi.find_dual_exhaustively().expect("dual should be found");
        let summary =
            summarize_dual_verification(&phi, &dual).expect("verification summary should build");

        assert!(summary.contains("duality kind: separable classical"));
        assert!(summary.contains("phi_hat(phi(Q)) = [3]Q on all Q in E(Fp): yes"));
        assert!(summary.contains("phi(phi_hat(R)) = [3]R on all R in E'(Fp): yes"));
        assert!(summary.contains("deg(phi_hat o phi) = 9 = deg([3]): yes"));
    }

    #[test]
    fn dual_report_description_mentions_kind_and_kernel_summaries() {
        let phi = degree_three_phi();
        let dual = phi.find_dual_exhaustively().expect("dual should be found");
        let report =
            crate::elliptic_curves::short_weierstrass::isogenies::DualVeluIsogeny::dual_report(
                &phi, &dual,
            )
            .expect("report should build");
        let description = describe_dual_isogeny_report(&report);

        assert!(description.contains("Dual isogeny report"));
        assert!(description.contains("duality kind: separable classical"));
        assert!(description.contains("phi kernel:"));
        assert!(description.contains("phi_hat kernel:"));
    }
}
