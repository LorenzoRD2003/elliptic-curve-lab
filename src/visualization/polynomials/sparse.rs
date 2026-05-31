use crate::fields::Field;
use crate::polynomials::{SparsePolynomial, SparsePolynomialTerm};
use crate::visualization::{Visualizable, VisualizableField, VisualizablePolynomial};

/// Formats a sparse polynomial as a human-readable univariate expression.
pub fn format_sparse_polynomial<F>(polynomial: &SparsePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut terms = Vec::new();

    for term in polynomial.terms().iter().rev() {
        let coefficient_text = term.coefficient.format_elem();
        let piece = match term.degree {
            0 => coefficient_text,
            1 if F::eq(&term.coefficient, &F::one()) => "x".to_string(),
            1 => format!("{coefficient_text}*x"),
            _ if F::eq(&term.coefficient, &F::one()) => format!("x^{}", term.degree),
            _ => format!("{coefficient_text}*x^{}", term.degree),
        };
        terms.push(piece);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

/// Explains the normalized sparse term list.
pub fn explain_sparse_storage<F>(polynomial: &SparsePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Sparse polynomial".to_string(),
        "terms are stored in ascending degree order".to_string(),
        format!("polynomial: {}", format_sparse_polynomial(polynomial)),
        "stored terms:".to_string(),
    ];

    if polynomial.is_empty() {
        lines.push("- no stored terms: this is the zero polynomial".to_string());
        return lines.join("\n");
    }

    for SparsePolynomialTerm {
        coefficient,
        degree,
    } in polynomial.terms()
    {
        lines.push(format!(
            "- degree {degree}: coefficient {}",
            coefficient.format_elem()
        ));
    }

    lines.join("\n")
}

/// Returns a richer textual description of a sparse polynomial.
pub fn describe_sparse_polynomial<F>(polynomial: &SparsePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let terms = polynomial
        .terms()
        .iter()
        .map(|term| (term.degree, term.coefficient.format_elem()))
        .collect::<Vec<_>>();

    format!(
        "Sparse polynomial\n\
         stored term count: {}\n\
         degree: {:?}\n\
         normalized terms (ascending): {:?}\n\
         expression: {}",
        polynomial.len(),
        polynomial.degree(),
        terms,
        format_sparse_polynomial(polynomial)
    )
}

impl<F> Visualizable for SparsePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_sparse_polynomial(self)
    }

    fn describe(&self) -> String {
        describe_sparse_polynomial(self)
    }
}

impl<F> VisualizablePolynomial for SparsePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_polynomial(&self) -> String {
        format_sparse_polynomial(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::{Field, Fp};
    use crate::polynomials::{SparsePolynomial, SparsePolynomialTerm};
    use crate::visualization::VisualizablePolynomial;

    use super::{describe_sparse_polynomial, explain_sparse_storage, format_sparse_polynomial};

    type F17 = Fp<17>;

    fn term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::elem_from_u64(coefficient),
            degree,
        }
    }

    #[test]
    fn sparse_formatter_handles_zero_and_multiple_terms() {
        assert_eq!(
            format_sparse_polynomial(&SparsePolynomial::<F17>::new(Vec::new())),
            "0"
        );

        let polynomial = SparsePolynomial::<F17>::new(vec![term(3, 0), term(1, 2), term(5, 1)]);
        assert_eq!(
            format_sparse_polynomial(&polynomial),
            "x^2 + 5 (mod 17)*x + 3 (mod 17)"
        );
    }

    #[test]
    fn sparse_storage_explanation_mentions_order_and_terms() {
        let polynomial = SparsePolynomial::<F17>::new(vec![term(3, 0), term(5, 2)]);
        let explanation = explain_sparse_storage(&polynomial);

        assert!(explanation.contains("ascending degree order"));
        assert!(explanation.contains("degree 0: coefficient 3 (mod 17)"));
        assert!(explanation.contains("degree 2: coefficient 5 (mod 17)"));
    }

    #[test]
    fn sparse_description_contains_key_metadata() {
        let polynomial = SparsePolynomial::<F17>::new(vec![term(3, 0), term(5, 2)]);
        let description = describe_sparse_polynomial(&polynomial);

        assert!(description.contains("Sparse polynomial"));
        assert!(description.contains("stored term count: 2"));
        assert!(description.contains("degree: Some(2)"));
        assert!(description.contains("(0, \"3 (mod 17)\")"));
    }

    #[test]
    fn sparse_polynomial_visualizable_trait_reuses_core_helpers() {
        let polynomial = SparsePolynomial::<F17>::new(vec![term(3, 0), term(5, 2)]);

        assert_eq!(
            polynomial.format_polynomial(),
            "5 (mod 17)*x^2 + 3 (mod 17)"
        );
        assert!(
            polynomial
                .describe_polynomial()
                .contains("Sparse polynomial")
        );
    }
}
