use crate::polynomials::MultivariatePolynomial;
use crate::polynomials::multivariate::{Monomial, MultivariateTerm};
use crate::visualization::VisualizableField;
use crate::visualization::polynomials::traits::VisualizablePolynomial;
use crate::visualization::traits::Visualizable;
use crate::visualization::*;

/// Formats a monomial using variable names `x_0`, `x_1`, ...
pub fn format_monomial(monomial: &Monomial) -> String {
    let mut factors = Vec::new();

    for (index, exponent) in monomial.exponents.iter().copied().enumerate() {
        match exponent {
            0 => {}
            1 => factors.push(format!("x_{index}")),
            _ => factors.push(format!("x_{index}^{exponent}")),
        }
    }

    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

/// Formats a multivariate polynomial as a human-readable expression.
pub fn format_multivariate_polynomial<F>(polynomial: &MultivariatePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut pieces = Vec::new();

    for term in polynomial.terms().iter().rev() {
        let coefficient_text = term.coefficient.format_elem();
        let monomial_text = format_monomial(&term.monomial);

        let piece = if monomial_text == "1" {
            coefficient_text
        } else if F::eq(&term.coefficient, &F::one()) {
            monomial_text
        } else {
            format!("{coefficient_text}*{monomial_text}")
        };

        pieces.push(piece);
    }

    if pieces.is_empty() {
        "0".to_string()
    } else {
        pieces.join(" + ")
    }
}

/// Explains the normalized multivariate term storage.
pub fn explain_multivariate_storage<F>(polynomial: &MultivariatePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Multivariate polynomial".to_string(),
        format!("arity: {}", polynomial.arity()),
        "terms are stored in ascending monomial order".to_string(),
        format!("polynomial: {}", format_multivariate_polynomial(polynomial)),
        "stored terms:".to_string(),
    ];

    if polynomial.is_empty() {
        lines.push("- no stored terms: this is the zero polynomial".to_string());
        return lines.join("\n");
    }

    for MultivariateTerm {
        coefficient,
        monomial,
    } in polynomial.terms()
    {
        lines.push(format!(
            "- monomial {}: coefficient {}",
            format_monomial(monomial),
            coefficient.format_elem()
        ));
    }

    lines.join("\n")
}

/// Returns a richer textual description of a multivariate polynomial.
pub fn describe_multivariate_polynomial<F>(polynomial: &MultivariatePolynomial<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    let terms = polynomial
        .terms()
        .iter()
        .map(|term| {
            (
                format_monomial(&term.monomial),
                term.coefficient.format_elem(),
            )
        })
        .collect::<Vec<_>>();

    format!(
        "Multivariate polynomial\n\
         arity: {}\n\
         stored term count: {}\n\
         max total degree: {:?}\n\
         normalized terms (ascending): {:?}\n\
         expression: {}",
        polynomial.arity(),
        polynomial.len(),
        polynomial.degree(),
        terms,
        format_multivariate_polynomial(polynomial)
    )
}

impl<F> Visualizable for MultivariatePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_compact(&self) -> String {
        format_multivariate_polynomial(self)
    }

    fn describe(&self) -> String {
        describe_multivariate_polynomial(self)
    }
}

impl<F> VisualizablePolynomial for MultivariatePolynomial<F>
where
    F: Field,
    F::Elem: VisualizableField,
{
    fn format_polynomial(&self) -> String {
        format_multivariate_polynomial(self)
    }
}

#[cfg(test)]
mod tests {

    use crate::polynomials::MultivariatePolynomial;
    use crate::polynomials::multivariate::{Monomial, MultivariateTerm};
    use crate::visualization::VisualizablePolynomial;

    use crate::visualization::polynomials::{
        describe_multivariate_polynomial, explain_multivariate_storage, format_monomial,
        format_multivariate_polynomial,
    };

    type F17 = crate::fields::Fp17;

    fn term(coefficient: u64, exponents: &[usize]) -> MultivariateTerm<F17> {
        MultivariateTerm {
            coefficient: F17::from_i64(coefficient),
            monomial: Monomial::new(exponents.to_vec()),
        }
    }

    #[test]
    fn monomial_formatter_handles_constants_and_mixed_exponents() {
        assert_eq!(format_monomial(&Monomial::new(vec![0, 0, 0])), "1");
        assert_eq!(format_monomial(&Monomial::new(vec![2, 0, 1])), "x_0^2*x_2");
    }

    #[test]
    fn multivariate_formatter_handles_zero_and_multiple_terms() {
        assert_eq!(
            format_multivariate_polynomial(
                &MultivariatePolynomial::<F17>::new(2, Vec::new()).unwrap()
            ),
            "0"
        );

        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![term(3, &[0, 0]), term(1, &[2, 0]), term(5, &[0, 1])],
        )
        .expect("polynomial should exist");

        assert_eq!(
            format_multivariate_polynomial(&polynomial),
            "x_0^2 + 5*x_1 + 3"
        );
    }

    #[test]
    fn multivariate_storage_explanation_mentions_arity_and_terms() {
        let polynomial =
            MultivariatePolynomial::<F17>::new(2, vec![term(3, &[0, 0]), term(5, &[1, 1])])
                .expect("polynomial should exist");
        let explanation = explain_multivariate_storage(&polynomial);

        assert!(explanation.contains("arity: 2"));
        assert!(explanation.contains("monomial 1: coefficient 3"));
        assert!(explanation.contains("monomial x_0*x_1: coefficient 5"));
    }

    #[test]
    fn multivariate_description_contains_key_metadata() {
        let polynomial =
            MultivariatePolynomial::<F17>::new(2, vec![term(3, &[0, 0]), term(5, &[1, 1])])
                .expect("polynomial should exist");
        let description = describe_multivariate_polynomial(&polynomial);

        assert!(description.contains("Multivariate polynomial"));
        assert!(description.contains("arity: 2"));
        assert!(description.contains("stored term count: 2"));
        assert!(description.contains("max total degree: Some(2)"));
        assert!(description.contains("(\"1\", \"3\")"));
    }

    #[test]
    fn multivariate_polynomial_visualizable_trait_reuses_core_helpers() {
        let polynomial =
            MultivariatePolynomial::<F17>::new(2, vec![term(3, &[0, 0]), term(5, &[1, 1])])
                .expect("polynomial should exist");

        assert_eq!(polynomial.format_polynomial(), "5*x_0*x_1 + 3");
        assert!(
            polynomial
                .describe_polynomial()
                .contains("Multivariate polynomial")
        );
    }
}
