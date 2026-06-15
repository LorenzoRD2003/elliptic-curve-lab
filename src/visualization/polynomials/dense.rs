use crate::fields::traits::Field;
use crate::polynomials::DensePolynomial;
use crate::visualization::fields::traits::VisualizableField;
use crate::visualization::polynomials::traits::VisualizablePolynomial;
use crate::visualization::traits::Visualizable;

/// Formats a dense polynomial as a human-readable univariate expression.
pub fn format_dense_polynomial<F>(polynomial: &DensePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut terms = Vec::new();

    for (degree, coefficient) in polynomial.coefficients().iter().enumerate().rev() {
        if F::is_zero(coefficient) {
            continue;
        }

        let coefficient_text = coefficient.format_elem();
        let term = match degree {
            0 => coefficient_text,
            1 if F::eq(coefficient, &F::one()) => "x".to_string(),
            1 => format!("{coefficient_text}*x"),
            _ if F::eq(coefficient, &F::one()) => format!("x^{degree}"),
            _ => format!("{coefficient_text}*x^{degree}"),
        };

        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

/// Explains how the dense coefficient vector maps to a polynomial.
pub fn explain_dense_storage<F>(polynomial: &DensePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Dense polynomial".to_string(),
        "coefficients are stored in ascending degree order".to_string(),
        format!("polynomial: {}", format_dense_polynomial(polynomial)),
        "storage mapping:".to_string(),
    ];

    if polynomial.is_zero() {
        lines.push("- empty vector represents the zero polynomial".to_string());
        return lines.join("\n");
    }

    for (degree, coefficient) in polynomial.coefficients().iter().enumerate() {
        lines.push(format!(
            "- index {degree}: coefficient {} multiplies x^{degree}",
            coefficient.format_elem()
        ));
    }

    lines.join("\n")
}

/// Returns a richer textual description of a dense polynomial.
pub fn describe_dense_polynomial<F>(polynomial: &DensePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let coefficients = polynomial
        .coefficients()
        .iter()
        .map(VisualizableField::format_elem)
        .collect::<Vec<_>>();

    format!(
        "Dense polynomial\n\
         coefficient count: {}\n\
         degree: {:?}\n\
         raw coefficients (ascending): {:?}\n\
         expression: {}",
        polynomial.len(),
        polynomial.degree(),
        coefficients,
        format_dense_polynomial(polynomial)
    )
}

impl<F> Visualizable for DensePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_dense_polynomial(self)
    }

    fn describe(&self) -> String {
        describe_dense_polynomial(self)
    }
}

impl<F> VisualizablePolynomial for DensePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_polynomial(&self) -> String {
        format_dense_polynomial(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::{Fp, traits::Field};
    use crate::polynomials::DensePolynomial;
    use crate::visualization::VisualizablePolynomial;

    use crate::visualization::polynomials::{
        describe_dense_polynomial, explain_dense_storage, format_dense_polynomial,
    };

    type F17 = Fp<17>;

    fn coeffs(values: &[u64]) -> Vec<<F17 as Field>::Elem> {
        values.iter().copied().map(F17::elem_from_u64).collect()
    }

    #[test]
    fn dense_formatter_handles_zero_and_sparse_internal_coefficients() {
        assert_eq!(
            format_dense_polynomial(&DensePolynomial::<F17>::new(coeffs(&[]))),
            "0"
        );
        assert_eq!(
            format_dense_polynomial(&DensePolynomial::<F17>::new(coeffs(&[3, 0, 1]))),
            "x^2 + 3"
        );
    }

    #[test]
    fn dense_storage_explanation_mentions_order_and_indices() {
        let polynomial = DensePolynomial::<F17>::new(coeffs(&[3, 0, 1]));
        let explanation = explain_dense_storage(&polynomial);

        assert!(explanation.contains("ascending degree order"));
        assert!(explanation.contains("index 0: coefficient 3 multiplies x^0"));
        assert!(explanation.contains("index 2: coefficient 1 multiplies x^2"));
    }

    #[test]
    fn dense_description_contains_key_metadata() {
        let polynomial = DensePolynomial::<F17>::new(coeffs(&[2, 5, 1]));
        let description = describe_dense_polynomial(&polynomial);

        assert!(description.contains("Dense polynomial"));
        assert!(description.contains("coefficient count: 3"));
        assert!(description.contains("degree: Some(2)"));
        assert!(description.contains("expression: x^2 + 5*x + 2"));
    }

    #[test]
    fn dense_polynomial_visualizable_trait_reuses_core_helpers() {
        let polynomial = DensePolynomial::<F17>::new(coeffs(&[2, 5, 1]));

        assert_eq!(polynomial.format_polynomial(), "x^2 + 5*x + 2");
        assert!(
            polynomial
                .describe_polynomial()
                .contains("Dense polynomial")
        );
    }
}
