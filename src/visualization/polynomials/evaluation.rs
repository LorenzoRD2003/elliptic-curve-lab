use num_bigint::BigUint;

use crate::fields::traits::Field;
use crate::polynomials::{
    DensePolynomial, MultivariatePolynomial, PolynomialError, SparsePolynomial,
};
use crate::visualization::{
    VisualizableField,
    polynomials::{
        dense::format_dense_polynomial,
        multivariate::{format_monomial, format_multivariate_polynomial},
        sparse::format_sparse_polynomial,
    },
    shared::comma_list,
};

/// Explains dense univariate evaluation using Horner's rule.
fn explain_evaluate_dense<F: Field>(polynomial: &DensePolynomial<F>, point: &F::Elem) -> String
where
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Dense polynomial evaluation".to_string(),
        format!("polynomial: {}", format_dense_polynomial(polynomial)),
        format!("point: {}", point.format_elem()),
        "method: Horner's rule over coefficients in descending degree order".to_string(),
    ];

    let mut accumulator = F::zero();

    if polynomial.is_zero() {
        lines.push("the polynomial is zero, so the result is 0".to_string());
    } else {
        for (step, coefficient) in polynomial.coefficients().iter().rev().enumerate() {
            let previous = accumulator.clone();
            let multiplied = F::mul(&previous, point);
            accumulator = F::add(&multiplied, coefficient);

            lines.push(format!(
                "step {step}: ({}) * {} + {} = {}",
                previous.format_elem(),
                point.format_elem(),
                coefficient.format_elem(),
                accumulator.format_elem()
            ));
        }
    }

    let value = polynomial
        .evaluate(point)
        .expect("dense evaluation should not fail");
    lines.push(format!("result: {}", value.format_elem()));
    lines.join("\n")
}

/// Explains sparse univariate evaluation by summing explicit term
/// contributions.
fn explain_evaluate_sparse<F: Field>(polynomial: &SparsePolynomial<F>, point: &F::Elem) -> String
where
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Sparse polynomial evaluation".to_string(),
        format!("polynomial: {}", format_sparse_polynomial(polynomial)),
        format!("point: {}", point.format_elem()),
        "method: evaluate each stored non-zero term and add the contributions".to_string(),
    ];

    let mut total = F::zero();

    if polynomial.is_empty() {
        lines.push("the polynomial has no stored terms, so the result is 0".to_string());
    } else {
        for term in polynomial.terms() {
            let power = F::pow(point, &BigUint::from(term.degree));
            let contribution = F::mul(&term.coefficient, &power);
            total = F::add(&total, &contribution);

            lines.push(format!(
                "degree {} term: {} * {}^{} = {}",
                term.degree,
                term.coefficient.format_elem(),
                point.format_elem(),
                term.degree,
                contribution.format_elem()
            ));
        }
    }

    let value = polynomial
        .evaluate(point)
        .expect("sparse evaluation should not fail");
    lines.push(format!("result: {}", value.format_elem()));
    lines.join("\n")
}

/// Explains multivariate evaluation by evaluating each monomial at the chosen
/// point and summing the term contributions.
fn explain_evaluate_multivariate<F: Field>(
    polynomial: &MultivariatePolynomial<F>,
    point: &[F::Elem],
) -> Result<String, PolynomialError>
where
    F::Elem: VisualizableField,
{
    if point.len() != polynomial.arity() {
        return Err(PolynomialError::EvaluationPointArityMismatch {
            expected: polynomial.arity(),
            actual: point.len(),
        });
    }

    let mut lines = vec![
        "Multivariate polynomial evaluation".to_string(),
        format!("polynomial: {}", format_multivariate_polynomial(polynomial)),
        format!(
            "point: ({})",
            comma_list(point.iter().map(VisualizableField::format_elem))
        ),
        "method: evaluate each monomial coordinate-wise and sum the term contributions".to_string(),
    ];

    let mut total = F::zero();

    if polynomial.is_empty() {
        lines.push("the polynomial has no stored terms, so the result is 0".to_string());
    } else {
        for term in polynomial.terms() {
            let mut monomial_value = F::one();

            for (coordinate, exponent) in point.iter().zip(term.monomial().exponents()) {
                let power = F::pow(coordinate, &BigUint::from(*exponent));
                monomial_value = F::mul(&monomial_value, &power);
            }

            let contribution = F::mul(term.coefficient(), &monomial_value);
            total = F::add(&total, &contribution);

            lines.push(format!(
                "term {} with coefficient {} contributes {}",
                format_monomial(term.monomial()),
                term.coefficient().format_elem(),
                contribution.format_elem()
            ));
        }
    }

    let value = polynomial.evaluate(point)?;
    lines.push(format!("result: {}", value.format_elem()));
    Ok(lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomials::{
        DensePolynomial, MultivariatePolynomial, PolynomialError, SparsePolynomial,
        multivariate::{Monomial, MultivariateTerm},
        sparse::SparsePolynomialTerm,
    };

    type F17 = crate::fields::Fp17;

    fn dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::new(values.iter().copied().map(F17::from_i64).collect())
    }

    fn sparse_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::from_i64(coefficient),
            degree,
        }
    }

    fn multivariate_term(coefficient: u64, exponents: &[usize]) -> MultivariateTerm<F17> {
        MultivariateTerm::new(
            F17::from_i64(coefficient),
            Monomial::new(exponents.to_vec()),
        )
    }

    #[test]
    fn dense_evaluation_explanation_mentions_horner_and_result() {
        let explanation = explain_evaluate_dense(&dense(&[3, 5, 2]), &F17::from_i64(4));
        assert!(explanation.contains("Dense polynomial evaluation"));
        assert!(explanation.contains("Horner"));
        assert!(explanation.contains("result: 4"));
    }

    #[test]
    fn sparse_evaluation_explanation_mentions_term_contributions() {
        let polynomial = SparsePolynomial::<F17>::new(vec![
            sparse_term(3, 0),
            sparse_term(5, 2),
            sparse_term(1, 3),
        ]);
        let explanation = explain_evaluate_sparse(&polynomial, &F17::from_i64(2));
        assert!(explanation.contains("Sparse polynomial evaluation"));
        assert!(explanation.contains("degree 2 term"));
        assert!(explanation.contains("result: 14"));
    }

    #[test]
    fn multivariate_evaluation_explanation_mentions_monomials_and_result() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![
                multivariate_term(3, &[0, 0]),
                multivariate_term(5, &[1, 1]),
                multivariate_term(1, &[2, 0]),
            ],
        )
        .expect("polynomial should exist");
        let point = [F17::from_i64(2), F17::from_i64(3)];
        let explanation =
            explain_evaluate_multivariate(&polynomial, &point).expect("evaluation should work");

        assert!(explanation.contains("Multivariate polynomial evaluation"));
        assert!(explanation.contains("term x_0*x_1"));
        assert!(explanation.contains("result: 3"));
    }

    #[test]
    fn multivariate_evaluation_explanation_rejects_wrong_arity() {
        let polynomial = MultivariatePolynomial::<F17>::new(2, vec![multivariate_term(1, &[1, 0])])
            .expect("polynomial should exist");
        let point = [F17::from_i64(2)];

        let error = explain_evaluate_multivariate(&polynomial, &point)
            .expect_err("wrong arity should fail");

        assert_eq!(
            error,
            PolynomialError::EvaluationPointArityMismatch {
                expected: 2,
                actual: 1,
            }
        );
    }
}
