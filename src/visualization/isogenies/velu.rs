use core::fmt;
use std::hash::Hash;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{
    CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::fields::{EnumerableFiniteField, Field, FiniteField, SqrtField};
use crate::isogenies::{Isogeny, IsogenyError, VeluIsogeny};
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::traits::Visualizable;

use crate::visualization::elliptic_curves::{format_curve, format_point, format_point_compact};

fn format_elem<F>(value: &F::Elem) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    value.format_elem()
}

fn field_surface<F>() -> String
where
    F: FiniteField,
{
    let characteristic = F::characteristic();
    let extension_degree = F::extension_degree().get();
    if extension_degree == 1 {
        format!("F_{characteristic}")
    } else {
        format!("F_{}^{}", characteristic, extension_degree)
    }
}

fn compact_curve_surface<F>(curve: &ShortWeierstrassCurve<F>) -> String
where
    F: FiniteField,
    F::Elem: VisualizableField + fmt::Display,
{
    format!(
        "y^2 = x^3 + {}x + {} over {}",
        format_elem::<F>(curve.a()),
        format_elem::<F>(curve.b()),
        field_surface::<F>()
    )
}

fn format_kernel_points<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> String
where
    F: Field + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    isogeny
        .kernel_points()
        .iter()
        .map(format_point_compact::<F>)
        .collect::<Vec<_>>()
        .join(", ")
}

fn cyclic_kernel_generator<F>(
    isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>,
) -> Option<&AffinePoint<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    isogeny
        .kernel_nonzero_points()
        .iter()
        .find(|point| isogeny.domain().point_order(point) == Some(isogeny.kernel().order()))
}

fn kernel_maps_to_identity<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> bool
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    let identity = isogeny.codomain().identity();
    isogeny.kernel_points().iter().all(|point| {
        isogeny
            .evaluate(point)
            .ok()
            .is_some_and(|image| image == identity)
    })
}

fn constant_on_kernel_cosets<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> bool
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    let domain_points = isogeny.domain().points();
    domain_points.iter().all(|point| {
        let image = match isogeny.evaluate(point) {
            Ok(image) => image,
            Err(_) => return false,
        };

        isogeny.kernel_points().iter().all(|kernel_point| {
            let translated = match isogeny.domain().add(point, kernel_point) {
                Ok(translated) => translated,
                Err(_) => return false,
            };

            isogeny
                .evaluate(&translated)
                .ok()
                .is_some_and(|translated_image| translated_image == image)
        })
    })
}

fn exhaustive_homomorphism_check<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> bool
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    let domain_points = isogeny.domain().points();
    domain_points.iter().all(|left| {
        domain_points.iter().all(|right| {
            let domain_sum = match isogeny.domain().add(left, right) {
                Ok(sum) => sum,
                Err(_) => return false,
            };
            let image_left = match isogeny.evaluate(left) {
                Ok(image) => image,
                Err(_) => return false,
            };
            let image_right = match isogeny.evaluate(right) {
                Ok(image) => image,
                Err(_) => return false,
            };
            let image_sum = match isogeny.evaluate(&domain_sum) {
                Ok(image) => image,
                Err(_) => return false,
            };
            let codomain_sum = match isogeny.codomain().add(&image_left, &image_right) {
                Ok(sum) => sum,
                Err(_) => return false,
            };

            image_sum == codomain_sum
        })
    })
}

/// Formats a short-Weierstrass Vélu isogeny compactly.
pub fn format_isogeny<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> String
where
    F: FiniteField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    format!(
        "Velu isogeny of degree {}: {} -> {}",
        isogeny.degree(),
        compact_curve_surface(isogeny.domain()),
        compact_curve_surface(isogeny.codomain())
    )
}

/// Describes the main structural data of a short-Weierstrass Vélu isogeny.
pub fn describe_isogeny<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> String
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    let kernel_header = match cyclic_kernel_generator(isogeny) {
        Some(generator) => format!(
            "  order {} cyclic subgroup generated by P = {}",
            isogeny.kernel().order(),
            format_point_compact(generator)
        ),
        None => format!("  finite subgroup of order {}", isogeny.kernel().order()),
    };

    let mut lines = vec![
        "Vélu isogeny".to_string(),
        "============".to_string(),
        String::new(),
        "domain:".to_string(),
        format!("  E: {}", compact_curve_surface(isogeny.domain())),
        String::new(),
        "kernel:".to_string(),
        kernel_header,
        "  points:".to_string(),
    ];

    lines.extend(
        isogeny
            .kernel_points()
            .iter()
            .map(|point| format!("    {}", format_point_compact(point))),
    );

    lines.extend([
        String::new(),
        "codomain:".to_string(),
        format!("  E': {}", compact_curve_surface(isogeny.codomain())),
        String::new(),
        "degree:".to_string(),
        format!("  {}", isogeny.degree()),
        String::new(),
        "checks:".to_string(),
        format!(
            "  kernel maps to identity: {}",
            if kernel_maps_to_identity(isogeny) {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "  φ(P + Q) = φ(P) for Q in kernel: {}",
            if constant_on_kernel_cosets(isogeny) {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "  homomorphism check on all points: {}",
            if exhaustive_homomorphism_check(isogeny) {
                "yes"
            } else {
                "no"
            }
        ),
        format!(
            "  #E({}) = #E'({}): {}",
            field_surface::<F>(),
            field_surface::<F>(),
            if isogeny.domain().order() == isogeny.codomain().order() {
                "yes"
            } else {
                "no"
            }
        ),
    ]);

    lines.join("\n")
}

/// Summarizes the explicit kernel subgroup carried by a Vélu isogeny.
pub fn summarize_kernel<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> String
where
    F: FiniteField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    [
        "Kernel summary".to_string(),
        format!("order: {}", isogeny.kernel().order()),
        format!("degree contribution: {}", isogeny.kernel().degree()),
        format!(
            "non-identity points: {}",
            isogeny.kernel_nonzero_points().len()
        ),
        format!("points: {}", format_kernel_points(isogeny)),
    ]
    .join("\n")
}

/// Explains the codomain coefficient computation for a short-Weierstrass Vélu isogeny.
pub fn explain_velu_codomain<F>(isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>) -> String
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    let domain = isogeny.domain();
    let a = domain.a().clone();
    let b = domain.b().clone();
    let mut a_sum = F::zero();
    let mut b_sum = F::zero();
    let mut lines = vec![
        "Vélu codomain explanation".to_string(),
        format!("domain: {}", format_curve(domain)),
        format!("kernel: {}", format_kernel_points(isogeny)),
        "formula: a' = a - 5 * sum_{Q in G*} (3x_Q^2 + a)".to_string(),
        "formula: b' = b - 7 * sum_{Q in G*} (5x_Q^3 + 3ax_Q + 2b)".to_string(),
    ];

    for kernel_point in isogeny.kernel_nonzero_points() {
        let x = AffinePoint::x_coordinate(kernel_point)
            .expect("non-identity kernel points should be finite in the current affine model");
        let x_squared = F::square(x);
        let x_cubed = F::cube(x);
        let a_term = F::add(&F::mul(&F::from_i64(3), &x_squared), &a);
        let b_term = F::add(
            &F::add(
                &F::mul(&F::from_i64(5), &x_cubed),
                &F::mul(&F::from_i64(3), &F::mul(&a, x)),
            ),
            &F::mul(&F::from_i64(2), &b),
        );

        a_sum = F::add(&a_sum, &a_term);
        b_sum = F::add(&b_sum, &b_term);

        lines.push(format!("Q: {}", format_point(kernel_point)));
        lines.push(format!("x_Q: {}", format_elem::<F>(x)));
        lines.push(format!(
            "a-term: 3x_Q^2 + a = {} + {} = {}",
            format_elem::<F>(&F::mul(&F::from_i64(3), &x_squared)),
            format_elem::<F>(&a),
            format_elem::<F>(&a_term)
        ));
        lines.push(format!(
            "b-term: 5x_Q^3 + 3ax_Q + 2b = {} + {} + {} = {}",
            format_elem::<F>(&F::mul(&F::from_i64(5), &x_cubed)),
            format_elem::<F>(&F::mul(&F::from_i64(3), &F::mul(&a, x))),
            format_elem::<F>(&F::mul(&F::from_i64(2), &b)),
            format_elem::<F>(&b_term)
        ));
    }

    let a_prime = F::sub(&a, &F::mul(&F::from_i64(5), &a_sum));
    let b_prime = F::sub(&b, &F::mul(&F::from_i64(7), &b_sum));

    lines.push(format!("sum_a: {}", format_elem::<F>(&a_sum)));
    lines.push(format!("sum_b: {}", format_elem::<F>(&b_sum)));
    lines.push(format!(
        "a': {} - 5 * {} = {}",
        format_elem::<F>(&a),
        format_elem::<F>(&a_sum),
        format_elem::<F>(&a_prime)
    ));
    lines.push(format!(
        "b': {} - 7 * {} = {}",
        format_elem::<F>(&b),
        format_elem::<F>(&b_sum),
        format_elem::<F>(&b_prime)
    ));
    lines.push(format!("codomain: {}", format_curve(isogeny.codomain())));

    lines.join("\n")
}

/// Explains the current short-Weierstrass Vélu evaluation on one point.
pub fn explain_velu_evaluation<F>(
    isogeny: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    point: &AffinePoint<F>,
) -> Result<String, IsogenyError>
where
    F: Field + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    if !isogeny.domain().contains(point) {
        return Err(IsogenyError::Curve(CurveError::PointNotOnCurve));
    }

    let mut lines = vec![
        "Vélu evaluation explanation".to_string(),
        format!("domain: {}", format_curve(isogeny.domain())),
        format!("codomain: {}", format_curve(isogeny.codomain())),
        format!("point: {}", format_point(point)),
    ];

    if isogeny.kernel().contains(point) {
        lines.push("case: the point lies in the kernel".to_string());
        lines.push("result: phi(P) = O".to_string());
        return Ok(lines.join("\n"));
    }

    let Some((x_phi, y_phi)) = isogeny.translation_sum_coordinates(point)? else {
        lines.push(
            "case: the image is the identity, so there are no affine output coordinates"
                .to_string(),
        );
        lines.push("result: phi(P) = O".to_string());
        return Ok(lines.join("\n"));
    };

    lines.push("formula: x_phi(P) = x(P) + sum_{Q in G*} (x(P + Q) - x(Q))".to_string());
    lines.push("formula: y_phi(P) = y(P) + sum_{Q in G*} (y(P + Q) - y(Q))".to_string());

    let x = AffinePoint::x_coordinate(point)
        .expect("a non-kernel affine point should have an x-coordinate");
    let y = AffinePoint::y_coordinate(point)
        .expect("a non-kernel affine point should have a y-coordinate");
    let mut x_sum = x.clone();
    let mut y_sum = y.clone();

    lines.push(format!("initial x(P): {}", format_elem::<F>(x)));
    lines.push(format!("initial y(P): {}", format_elem::<F>(y)));

    for kernel_point in isogeny.kernel_nonzero_points() {
        let translated = isogeny.domain().add(point, kernel_point)?;
        let (translated_x, translated_y) = AffinePoint::finite_coordinates(&translated)
            .expect("a non-kernel point translated by a non-zero kernel point should stay finite");
        let (kernel_x, kernel_y) = AffinePoint::finite_coordinates(kernel_point)
            .expect("non-identity kernel points should be finite");
        let x_delta = F::sub(translated_x, kernel_x);
        let y_delta = F::sub(translated_y, kernel_y);
        x_sum = F::add(&x_sum, &x_delta);
        y_sum = F::add(&y_sum, &y_delta);

        lines.push(format!("Q: {}", format_point(kernel_point)));
        lines.push(format!("P + Q: {}", format_point(&translated)));
        lines.push(format!(
            "x(P + Q) - x(Q): {} - {} = {}",
            format_elem::<F>(translated_x),
            format_elem::<F>(kernel_x),
            format_elem::<F>(&x_delta)
        ));
        lines.push(format!(
            "y(P + Q) - y(Q): {} - {} = {}",
            format_elem::<F>(translated_y),
            format_elem::<F>(kernel_y),
            format_elem::<F>(&y_delta)
        ));
    }

    let image = isogeny.evaluate(point)?;
    lines.push(format!("x_phi(P): {}", format_elem::<F>(&x_phi)));
    lines.push(format!("y_phi(P): {}", format_elem::<F>(&y_phi)));
    lines.push(format!("affine image: {}", format_point(&image)));
    debug_assert!(isogeny.codomain().contains(&image));
    Ok(lines.join("\n"))
}

impl<F> Visualizable for VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: VisualizableField + fmt::Display + Clone + Eq + Hash,
{
    fn format_compact(&self) -> String {
        format_isogeny(self)
    }

    fn describe(&self) -> String {
        describe_isogeny(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        AffineCurveModel, AffinePoint, CurveError, ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::IsogenyError;
    use crate::visualization::Visualizable;

    use crate::isogenies::VeluIsogeny;
    use crate::visualization::isogenies::{
        describe_isogeny, explain_velu_codomain, explain_velu_evaluation, format_isogeny,
        summarize_kernel,
    };

    type F41 = Fp<41>;

    fn f41_curve() -> ShortWeierstrassCurve<F41> {
        ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn f41_isogeny() -> VeluIsogeny<ShortWeierstrassCurve<F41>> {
        let domain = f41_curve();
        let generator = domain
            .point(F41::from_i64(40), F41::from_i64(0))
            .expect("point should lie on the curve");
        VeluIsogeny::from_generator(domain, generator).expect("Vélu isogeny should build")
    }

    #[test]
    fn format_isogeny_reports_degree_and_curves() {
        let formatted = format_isogeny(&f41_isogeny());

        assert!(formatted.contains("Velu isogeny of degree 2"));
        assert!(formatted.contains("->"));
        assert!(formatted.contains("y^2 = x^3 + 2x + 3 over F_41"));
        assert!(formatted.contains("y^2 = x^3 + 18x + 38 over F_41"));
    }

    #[test]
    fn describe_isogeny_mentions_kernel_and_construction() {
        let description = describe_isogeny(&f41_isogeny());

        assert!(description.contains("Vélu isogeny"));
        assert!(description.contains("domain:"));
        assert!(description.contains("E: y^2 = x^3 + 2x + 3 over F_41"));
        assert!(description.contains("order 2 cyclic subgroup generated by P = (40, 0)"));
        assert!(description.contains("    O"));
        assert!(description.contains("    (40, 0)"));
        assert!(description.contains("checks:"));
    }

    #[test]
    fn summarize_kernel_lists_the_small_explicit_subgroup() {
        let summary = summarize_kernel(&f41_isogeny());

        assert!(summary.contains("Kernel summary"));
        assert!(summary.contains("order: 2"));
        assert!(summary.contains("degree contribution: 2"));
        assert!(summary.contains("points: O, (40, 0)"));
    }

    #[test]
    fn velu_codomain_explanation_shows_the_updated_b_formula_and_result() {
        let explanation = explain_velu_codomain(&f41_isogeny());

        assert!(explanation.contains("Vélu codomain explanation"));
        assert!(explanation.contains("formula: a' = a - 5 * sum_{Q in G*} (3x_Q^2 + a)"));
        assert!(explanation.contains("formula: b' = b - 7 * sum_{Q in G*} (5x_Q^3 + 3ax_Q + 2b)"));
        assert!(explanation.contains("Q: (40 (mod 41), 0 (mod 41))"));
        assert!(explanation.contains("sum_a"));
        assert!(explanation.contains("sum_b"));
        assert!(explanation.contains("codomain: y^2 = x^3 + (18 (mod 41))x + (38 (mod 41))"));
    }

    #[test]
    fn velu_evaluation_explanation_reports_identity_for_kernel_points() {
        let isogeny = f41_isogeny();
        let explanation = explain_velu_evaluation(&isogeny, &AffinePoint::infinity())
            .expect("kernel identity explanation should succeed");

        assert!(explanation.contains("Vélu evaluation explanation"));
        assert!(explanation.contains("case: the point lies in the kernel"));
        assert!(explanation.contains("result: phi(P) = O"));
    }

    #[test]
    fn velu_evaluation_explanation_shows_translation_terms_for_non_kernel_points() {
        let isogeny = f41_isogeny();
        let point = f41_curve()
            .point(F41::from_i64(3), F41::from_i64(6))
            .expect("point should lie on the curve");
        let explanation =
            explain_velu_evaluation(&isogeny, &point).expect("evaluation explanation should work");

        assert!(explanation.contains("formula: x_phi(P) = x(P) + sum_{Q in G*} (x(P + Q) - x(Q))"));
        assert!(explanation.contains("formula: y_phi(P) = y(P) + sum_{Q in G*} (y(P + Q) - y(Q))"));
        assert!(explanation.contains("Q: (40 (mod 41), 0 (mod 41))"));
        assert!(explanation.contains("P + Q:"));
        assert!(explanation.contains("x(P + Q) - x(Q):"));
        assert!(explanation.contains("y(P + Q) - y(Q):"));
        assert!(explanation.contains("x_phi(P): 35"));
        assert!(explanation.contains("y_phi(P): 40"));
        assert!(explanation.contains("affine image: (35 (mod 41), 40 (mod 41))"));
    }

    #[test]
    fn velu_evaluation_explanation_rejects_off_curve_points() {
        let isogeny = f41_isogeny();
        let off_curve = AffinePoint::<F41>::new(F41::from_i64(2), F41::from_i64(2));

        assert_eq!(
            explain_velu_evaluation(&isogeny, &off_curve),
            Err(IsogenyError::Curve(CurveError::PointNotOnCurve))
        );
    }

    #[test]
    fn visualizable_trait_is_hooked_up_for_velu_isogenies() {
        let isogeny = f41_isogeny();

        assert!(
            isogeny
                .format_compact()
                .contains("Velu isogeny of degree 2")
        );
        assert!(isogeny.describe().contains("Vélu isogeny"));
    }
}
