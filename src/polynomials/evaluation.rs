use crate::fields::Field;
use crate::polynomials::{
    DensePolynomial, MultivariatePolynomial, PolynomialError, SparsePolynomial,
};

/// Evaluates a dense univariate polynomial over a field at a point of the same
/// field.
///
/// The implementation uses Horner's rule on the coefficient vector stored in
/// ascending degree order. This gives a clear and efficient baseline while
/// staying easy to explain.
///
/// If the polynomial is represented by an empty coefficient vector, the result
/// is the additive identity of the field.
pub fn evaluate_dense<F: Field>(
    polynomial: &DensePolynomial<F>,
    point: &F::Elem,
) -> Result<F::Elem, PolynomialError> {
    let mut accumulator = F::zero();

    for coefficient in polynomial.coefficients().iter().rev() {
        accumulator = F::add(&F::mul(&accumulator, point), coefficient);
    }

    Ok(accumulator)
}

/// Evaluates a sparse univariate polynomial over a field at a point of the
/// same field.
///
/// Because the representation stores only non-zero terms, the implementation
/// evaluates each stored term as
///
/// `coefficient * point^degree`
///
/// and sums the results.
pub fn evaluate_sparse<F: Field>(
    polynomial: &SparsePolynomial<F>,
    point: &F::Elem,
) -> Result<F::Elem, PolynomialError> {
    let mut value = F::zero();

    for term in polynomial.terms() {
        let power = F::pow(point, term.degree as u64);
        let contribution = F::mul(&term.coefficient, &power);
        value = F::add(&value, &contribution);
    }

    Ok(value)
}

/// Evaluates a multivariate polynomial over a field at a concrete point.
///
/// The point must provide exactly one coordinate for each variable in the
/// polynomial's ambient arity.
///
/// Each term is evaluated as
///
/// `coefficient * x_0^e0 * x_1^e1 * ... * x_{n-1}^{e_{n-1}}`
///
/// and the term contributions are then added together.
pub fn evaluate_multivariate<F: Field>(
    polynomial: &MultivariatePolynomial<F>,
    point: &[F::Elem],
) -> Result<F::Elem, PolynomialError> {
    if point.len() != polynomial.arity() {
        return Err(PolynomialError::EvaluationPointArityMismatch {
            expected: polynomial.arity(),
            actual: point.len(),
        });
    }

    let mut value = F::zero();

    for term in polynomial.terms() {
        let mut monomial_value = F::one();

        for (coordinate, exponent) in point.iter().zip(&term.monomial.exponents) {
            let power = F::pow(coordinate, *exponent as u64);
            monomial_value = F::mul(&monomial_value, &power);
        }

        let contribution = F::mul(&term.coefficient, &monomial_value);
        value = F::add(&value, &contribution);
    }

    Ok(value)
}

#[cfg(test)]
mod tests {
    use crate::fields::{Field, Fp, Q};
    use crate::polynomials::{
        DensePolynomial, Monomial, MultivariatePolynomial, MultivariateTerm, PolynomialError,
        SparsePolynomial, SparsePolynomialTerm,
    };

    use super::{evaluate_dense, evaluate_multivariate, evaluate_sparse};

    type F17 = Fp<17>;

    fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    fn f17_sparse_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::elem_from_u64(coefficient),
            degree,
        }
    }

    fn q(numerator: i64, denominator: i64) -> <Q as Field>::Elem {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        Q::div(&numerator, &denominator).expect("denominator should be non-zero")
    }

    fn q_sparse_term(numerator: i64, denominator: i64, degree: usize) -> SparsePolynomialTerm<Q> {
        SparsePolynomialTerm {
            coefficient: q(numerator, denominator),
            degree,
        }
    }

    fn q_multivariate_term(
        numerator: i64,
        denominator: i64,
        exponents: &[usize],
    ) -> MultivariateTerm<Q> {
        MultivariateTerm {
            coefficient: q(numerator, denominator),
            monomial: Monomial::new(exponents.to_vec()),
        }
    }

    #[test]
    fn evaluate_dense_uses_horner_rule_over_f17() {
        let polynomial = f17_dense(&[3, 5, 2]);
        let value =
            evaluate_dense(&polynomial, &F17::elem_from_u64(4)).expect("evaluation should work");

        assert!(F17::eq(&value, &F17::elem_from_u64(4)));
    }

    #[test]
    fn evaluate_dense_zero_polynomial_returns_zero() {
        let polynomial = DensePolynomial::<F17>::new(Vec::new());
        let value =
            evaluate_dense(&polynomial, &F17::elem_from_u64(9)).expect("evaluation should work");

        assert!(F17::eq(&value, &F17::zero()));
    }

    #[test]
    fn evaluate_sparse_matches_dense_over_f17() {
        let dense = f17_dense(&[3, 0, 5, 1]);
        let sparse = SparsePolynomial::<F17>::new(vec![
            f17_sparse_term(3, 0),
            f17_sparse_term(5, 2),
            f17_sparse_term(1, 3),
        ]);
        let point = F17::elem_from_u64(2);

        let dense_value = evaluate_dense(&dense, &point).expect("dense evaluation should work");
        let sparse_value = evaluate_sparse(&sparse, &point).expect("sparse evaluation should work");

        assert!(F17::eq(&dense_value, &sparse_value));
        assert!(F17::eq(&dense_value, &F17::elem_from_u64(14)));
    }

    #[test]
    fn evaluate_sparse_works_over_q() {
        let polynomial =
            SparsePolynomial::<Q>::new(vec![q_sparse_term(1, 2, 0), q_sparse_term(2, 3, 2)]);
        let value = evaluate_sparse(&polynomial, &q(3, 2)).expect("evaluation should work");

        assert!(Q::eq(&value, &q(2, 1)));
    }

    #[test]
    fn evaluate_multivariate_works_over_f17() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![
                MultivariateTerm {
                    coefficient: F17::elem_from_u64(3),
                    monomial: Monomial::new(vec![0, 0]),
                },
                MultivariateTerm {
                    coefficient: F17::elem_from_u64(5),
                    monomial: Monomial::new(vec![1, 1]),
                },
                MultivariateTerm {
                    coefficient: F17::elem_from_u64(1),
                    monomial: Monomial::new(vec![2, 0]),
                },
            ],
        )
        .expect("polynomial should exist");

        let point = [F17::elem_from_u64(2), F17::elem_from_u64(3)];
        let value = evaluate_multivariate(&polynomial, &point).expect("evaluation should work");

        assert!(F17::eq(&value, &F17::elem_from_u64(3)));
    }

    #[test]
    fn evaluate_multivariate_rejects_wrong_arity() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![MultivariateTerm {
                coefficient: F17::elem_from_u64(1),
                monomial: Monomial::new(vec![1, 0]),
            }],
        )
        .expect("polynomial should exist");

        let point = [F17::elem_from_u64(2)];
        let error = evaluate_multivariate(&polynomial, &point).expect_err("arity should fail");

        assert_eq!(
            error,
            PolynomialError::EvaluationPointArityMismatch {
                expected: 2,
                actual: 1,
            }
        );
    }

    #[test]
    fn evaluate_multivariate_works_over_q() {
        let polynomial = MultivariatePolynomial::<Q>::new(
            2,
            vec![
                q_multivariate_term(1, 2, &[1, 0]),
                q_multivariate_term(2, 3, &[0, 1]),
                q_multivariate_term(1, 4, &[0, 0]),
            ],
        )
        .expect("polynomial should exist");

        let point = [q(3, 2), q(6, 5)];
        let value = evaluate_multivariate(&polynomial, &point).expect("evaluation should work");

        assert!(Q::eq(&value, &q(9, 5)));
    }
}
