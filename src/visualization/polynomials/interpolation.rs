use crate::polynomials::DensePolynomial;
use crate::polynomials::PolynomialError;
use crate::visualization::VisualizableField;
use crate::visualization::*;

use crate::visualization::polynomials::format_dense_polynomial;

/// Explains classical Lagrange interpolation over a field.
///
/// The explanation is intentionally structural rather than fully expanding each
/// basis polynomial coefficient-by-coefficient. It focuses on:
///
/// - the input samples
/// - the existence of one basis polynomial per sample
/// - the final interpolated polynomial
/// - a verification pass showing that the result matches every sample
pub fn explain_lagrange_interpolation<F>(
    samples: &[(F::Elem, F::Elem)],
) -> Result<String, PolynomialError>
where
    F: Field,
    F::Elem: VisualizableField,
{
    let polynomial = DensePolynomial::<F>::lagrange_interpolate(samples)?;

    let mut lines = vec![
        "Lagrange interpolation".to_string(),
        format!("sample count: {}", samples.len()),
    ];

    if samples.is_empty() {
        lines.push("no samples: the interpolated polynomial is the zero polynomial".to_string());
        lines.push(format!("result: {}", format_dense_polynomial(&polynomial)));
        return Ok(lines.join("\n"));
    }

    lines.push("samples:".to_string());
    for (index, (x, y)) in samples.iter().enumerate() {
        lines.push(format!(
            "- sample {index}: ({}, {})",
            x.format_elem(),
            y.format_elem()
        ));
    }

    lines.push("basis idea:".to_string());
    lines.push(
        "for each sample i, construct a basis polynomial L_i(x) that is 1 at x_i and 0 at every other sample abscissa"
            .to_string(),
    );
    lines.push("the final interpolant is the sum of y_i * L_i(x) over all samples".to_string());
    lines.push(format!(
        "resulting polynomial: {}",
        format_dense_polynomial(&polynomial)
    ));

    lines.push("verification on the input samples:".to_string());
    for (index, (x, y)) in samples.iter().enumerate() {
        let value = polynomial.evaluate(x)?;
        lines.push(format!(
            "- p(x_{index}) = p({}) = {} (expected {})",
            x.format_elem(),
            value.format_elem(),
            y.format_elem()
        ));
    }

    Ok(lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use crate::fields::traits::*;

    use crate::fields::Q;
    use crate::polynomials::PolynomialError;

    use crate::visualization::polynomials::explain_lagrange_interpolation;

    type F17 = crate::fields::Fp17;

    fn q(numerator: i64, denominator: i64) -> <Q as crate::fields::traits::Field>::Elem {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        Q::div(&numerator, &denominator).expect("denominator should be non-zero")
    }

    #[test]
    fn lagrange_explanation_mentions_samples_and_result_over_f17() {
        let samples = [
            (F17::from_i64(0), F17::from_i64(3)),
            (F17::from_i64(1), F17::from_i64(10)),
            (F17::from_i64(2), F17::from_i64(4)),
        ];

        let explanation =
            explain_lagrange_interpolation::<F17>(&samples).expect("interpolation should work");

        assert!(explanation.contains("Lagrange interpolation"));
        assert!(explanation.contains("sample count: 3"));
        assert!(explanation.contains("sample 0: (0, 3)"));
        assert!(explanation.contains("resulting polynomial: 2*x^2 + 5*x + 3"));
        assert!(explanation.contains("verification on the input samples"));
    }

    #[test]
    fn lagrange_explanation_handles_empty_input() {
        let explanation =
            explain_lagrange_interpolation::<F17>(&[]).expect("empty interpolation should work");

        assert!(explanation.contains("no samples"));
        assert!(explanation.contains("result: 0"));
    }

    #[test]
    fn lagrange_explanation_works_over_q_too() {
        let samples = [(q(0, 1), q(1, 2)), (q(1, 1), q(7, 6)), (q(2, 1), q(17, 6))];

        let explanation =
            explain_lagrange_interpolation::<Q>(&samples).expect("interpolation should work");

        assert!(explanation.contains("sample 1: (1, 7/6)"));
        assert!(explanation.contains("resulting polynomial: 1/2*x^2 + 1/6*x + 1/2"));
        assert!(explanation.contains("expected 17/6"));
    }

    #[test]
    fn lagrange_explanation_rejects_duplicate_x_coordinates() {
        let samples = [
            (F17::from_i64(3), F17::from_i64(1)),
            (F17::from_i64(3), F17::from_i64(9)),
        ];

        let error = explain_lagrange_interpolation::<F17>(&samples)
            .expect_err("duplicate x values should fail");

        assert_eq!(error, PolynomialError::DuplicateInterpolationAbscissa);
    }
}
