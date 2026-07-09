use crate::visualization::*;
use core::fmt;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::visualization::{
    VisualizableField,
    elliptic_curves::short_weierstrass::{format_curve, format_point_compact},
    polynomials::dense::format_dense_polynomial,
    shared::{comma_list, format_field_elem as format_elem},
};

fn format_points<F: Field>(points: &[AffinePoint<F>]) -> String
where
    F::Elem: VisualizableField + fmt::Display,
{
    if points.is_empty() {
        "[]".to_string()
    } else {
        comma_list(points.iter().map(format_point_compact::<F>))
    }
}

fn format_xs<F: Field>(xs: &[F::Elem]) -> String
where
    F::Elem: VisualizableField,
{
    if xs.is_empty() {
        "[]".to_string()
    } else {
        comma_list(xs.iter().map(format_elem::<F>))
    }
}

/// Shape classification for a short-Weierstrass division polynomial.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DivisionPolynomialKind {
    InX,
    YTimesX,
}

impl DivisionPolynomialKind {
    fn as_text(self) -> &'static str {
        match self {
            Self::InX => "polinomio en x",
            Self::YTimesX => "y veces polinomio en x",
        }
    }
}

/// Compact division-polynomial summary for one division-polynomial index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DivisionPolynomialSummary {
    pub n: usize,
    pub form: DivisionPolynomialKind,
    pub degree_in_x: Option<usize>,
    pub rational_root_count: usize,
    pub rational_torsion_point_count: usize,
    pub exact_order_point_count: usize,
}

fn kind_of_form<F: Field>(form: &DivisionPolynomialForm<F>) -> DivisionPolynomialKind {
    match form {
        DivisionPolynomialForm::InX(_) => DivisionPolynomialKind::InX,
        DivisionPolynomialForm::YTimes(_) => DivisionPolynomialKind::YTimesX,
    }
}

fn expected_degree_in_x(n: usize) -> Option<usize> {
    match n {
        0 => None,
        1 => Some(0),
        _ if n.is_multiple_of(2) => Some((n * n).saturating_sub(4) / 2),
        _ => Some((n * n).saturating_sub(1) / 2),
    }
}

/// Summarizes the current division-polynomial and torsion picture for one
/// index `n`.
fn division_polynomial_summary<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<DivisionPolynomialSummary, DivisionPolynomialError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let form = curve.division_polynomial(n)?;
    let rational_roots = curve.rational_x_candidates_for_division_polynomial(n)?;
    let torsion_points = curve.torsion_points_from_division_polynomial(n)?;
    let exact_points = curve.exact_n_torsion_points_from_division_polynomial(n)?;

    Ok(DivisionPolynomialSummary {
        n,
        form: kind_of_form(&form),
        degree_in_x: expected_degree_in_x(n),
        rational_root_count: rational_roots.len(),
        rational_torsion_point_count: torsion_points.len(),
        exact_order_point_count: exact_points.len(),
    })
}

/// Explains the current division polynomial `ψ_n` for a short-Weierstrass
/// curve.
fn explain_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<String, DivisionPolynomialError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let form = curve.division_polynomial(n)?;
    let summary = division_polynomial_summary(curve, n)?;
    let polynomial_text = match &form {
        DivisionPolynomialForm::InX(polynomial) => format_dense_polynomial(polynomial),
        DivisionPolynomialForm::YTimes(polynomial) => {
            format!("y * ({})", format_dense_polynomial(polynomial))
        }
    };

    let lines = [
        "Division polynomial".to_string(),
        format!("Curva: {}", format_curve(curve)),
        format!("Índice: {}", n),
        format!("Forma de ψ_n: {}", summary.form.as_text()),
        format!(
            "Grado esperado en x: {}",
            summary
                .degree_in_x
                .map(|degree| degree.to_string())
                .unwrap_or_else(|| "indefinido para el polinomio cero".to_string())
        ),
        format!("Polinomio obtenido: {}", polynomial_text),
        format!(
            "Raíces racionales: {}",
            format_xs::<F>(&curve.rational_x_candidates_for_division_polynomial(n)?)
        ),
    ];

    Ok(lines.join("\n"))
}

/// Explains rational torsion recovery via the current division-polynomial
/// tooling.
fn explain_torsion_via_division_polynomial<F: EnumerableFiniteField + SqrtField>(
    curve: &ShortWeierstrassCurve<F>,
    n: usize,
) -> Result<String, DivisionPolynomialError>
where
    F::Elem: VisualizableField + fmt::Display,
{
    let form = curve.division_polynomial(n)?;
    let summary = division_polynomial_summary(curve, n)?;
    let rational_roots = curve.rational_x_candidates_for_division_polynomial(n)?;
    let lifted_candidates = curve.torsion_candidates_from_division_polynomial(n)?;
    let n_torsion_points = curve.torsion_points_from_division_polynomial(n)?;
    let exact_points = curve.exact_n_torsion_points_from_division_polynomial(n)?;
    let report = curve.compare_division_polynomial_torsion_with_enumeration(n)?;

    let polynomial_text = match &form {
        DivisionPolynomialForm::InX(polynomial) => format_dense_polynomial(polynomial),
        DivisionPolynomialForm::YTimes(polynomial) => {
            format!("y * ({})", format_dense_polynomial(polynomial))
        }
    };

    let lines = [
        "Torsion via division polynomial".to_string(),
        format!("Curva: {}", format_curve(curve)),
        format!("Índice: {}", n),
        format!("Forma de ψ_n: {}", summary.form.as_text()),
        format!(
            "Grado esperado en x: {}",
            summary
                .degree_in_x
                .map(|degree| degree.to_string())
                .unwrap_or_else(|| "indefinido para el polinomio cero".to_string())
        ),
        format!("Polinomio obtenido: {}", polynomial_text),
        format!("Raíces racionales: {}", format_xs::<F>(&rational_roots)),
        format!(
            "Puntos racionales levantados: {}",
            format_points::<F>(&lifted_candidates)
        ),
        format!(
            "Puntos que satisfacen [n]P = O: {}",
            format_points::<F>(&n_torsion_points)
        ),
        format!(
            "Puntos de orden exacto n: {}",
            format_points::<F>(&exact_points)
        ),
        "Comparación contra enumeración:".to_string(),
        format!(
            "  candidatos por polinomio: {}",
            report.polynomial_candidate_count()
        ),
        format!(
            "  puntos n-torsión por polinomio: {}",
            report.polynomial_n_torsion_count()
        ),
        format!(
            "  puntos n-torsión por enumeración: {}",
            report.enumerated_n_torsion_count()
        ),
        format!(
            "  puntos exactos por polinomio: {}",
            report.exact_order_polynomial_count()
        ),
        format!(
            "  puntos exactos por enumeración: {}",
            report.exact_order_enumerated_count()
        ),
        format!(
            "  faltantes desde el polinomio: {}",
            format_points::<F>(report.missing_from_polynomial())
        ),
        format!(
            "  extras desde el polinomio: {}",
            format_points::<F>(report.extra_from_polynomial())
        ),
    ];

    Ok(lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic_curves::ShortWeierstrassCurve;

    type F17 = crate::fields::Fp17;
    type F23 = crate::fields::Fp23;

    #[test]
    fn summary_reports_core_counts_for_order_three_example() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");
        let summary = division_polynomial_summary(&curve, 3).expect("summary should work");

        assert_eq!(summary.n, 3);
        assert_eq!(summary.form, DivisionPolynomialKind::InX);
        assert_eq!(summary.degree_in_x, Some(4));
        assert_eq!(summary.rational_root_count, 1);
        assert_eq!(summary.rational_torsion_point_count, 2);
        assert_eq!(summary.exact_order_point_count, 2);
    }

    #[test]
    fn summary_reports_even_shape_for_order_four_example() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");
        let summary = division_polynomial_summary(&curve, 4).expect("summary should work");

        assert_eq!(summary.form, DivisionPolynomialKind::YTimesX);
        assert_eq!(summary.degree_in_x, Some(6));
    }

    #[test]
    fn division_polynomial_explanation_mentions_requested_fields() {
        let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
            .expect("curve should be non-singular");
        let explanation = explain_division_polynomial(&curve, 3).expect("explanation should work");

        assert!(explanation.contains("Curva: y^2 = x^3"));
        assert!(explanation.contains("2"));
        assert!(explanation.contains("3"));
        assert!(explanation.contains("Índice: 3"));
        assert!(explanation.contains("Forma de ψ_n: polinomio en x"));
        assert!(explanation.contains("Grado esperado en x: 4"));
        assert!(explanation.contains("Polinomio obtenido:"));
        assert!(explanation.contains("Raíces racionales:"));
    }

    #[test]
    fn torsion_explanation_mentions_all_requested_sections() {
        let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
            .expect("curve should be non-singular");
        let explanation =
            explain_torsion_via_division_polynomial(&curve, 12).expect("explanation should work");

        assert!(explanation.contains("Curva: y^2 = x^3"));
        assert!(explanation.contains("2"));
        assert!(explanation.contains("3"));
        assert!(explanation.contains("Índice: 12"));
        assert!(explanation.contains("Forma de ψ_n: y veces polinomio en x"));
        assert!(explanation.contains("Grado esperado en x: 70"));
        assert!(explanation.contains("Polinomio obtenido:"));
        assert!(explanation.contains("Raíces racionales:"));
        assert!(explanation.contains("Puntos racionales levantados:"));
        assert!(explanation.contains("Puntos que satisfacen [n]P = O:"));
        assert!(explanation.contains("Puntos de orden exacto n:"));
        assert!(explanation.contains("Comparación contra enumeración:"));
    }
}
