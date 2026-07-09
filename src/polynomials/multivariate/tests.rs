#[cfg(test)]
#[allow(clippy::module_inception)]
mod tests {
    use crate::fields::traits::*;

    use proptest::prelude::*;

    use crate::fields::Q;
    use crate::polynomials::{
        MultivariatePolynomial, PolynomialError,
        multivariate::{Monomial, MultivariateTerm},
    };
    use crate::proptest_support::{
        config::PolynomialStrategyConfig, fields::arb_fp_elem,
        polynomials::arb_multivariate_polynomial,
    };

    type F17 = crate::fields::Fp17;

    fn f17_term(coefficient: u64, exponents: &[usize]) -> MultivariateTerm<F17> {
        MultivariateTerm::new(
            F17::from_i64(coefficient),
            Monomial::new(exponents.to_vec()),
        )
    }

    fn q_term(numerator: i64, denominator: i64, exponents: &[usize]) -> MultivariateTerm<Q> {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        MultivariateTerm::new(
            Q::div(&numerator, &denominator).expect("denominator should be non-zero"),
            Monomial::new(exponents.to_vec()),
        )
    }

    fn q(numerator: i64, denominator: i64) -> <Q as crate::fields::traits::Field>::Elem {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        Q::div(&numerator, &denominator).expect("denominator should be non-zero")
    }

    fn assert_f17_term(term: &MultivariateTerm<F17>, exponents: &[usize], coefficient: u64) {
        assert_eq!(term.monomial().exponents(), exponents);
        assert!(F17::eq(term.coefficient(), &F17::from_i64(coefficient)));
    }

    fn assert_q_term(
        term: &MultivariateTerm<Q>,
        exponents: &[usize],
        numerator: i64,
        denominator: i64,
    ) {
        assert_eq!(term.monomial().exponents(), exponents);
        assert!(Q::eq(term.coefficient(), &q(numerator, denominator)));
    }

    #[test]
    fn evaluate_multivariate_works_over_f17() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![
                MultivariateTerm::from_exponents(F17::from_i64(3), vec![0, 0]),
                MultivariateTerm::from_exponents(F17::from_i64(5), vec![1, 1]),
                MultivariateTerm::from_exponents(F17::from_i64(1), vec![2, 0]),
            ],
        )
        .expect("polynomial should exist");

        let point = [F17::from_i64(2), F17::from_i64(3)];
        let value = polynomial.evaluate(&point).expect("evaluation should work");

        assert!(F17::eq(&value, &F17::from_i64(3)));
    }

    #[test]
    fn evaluate_multivariate_rejects_wrong_arity() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![MultivariateTerm::from_exponents(
                F17::from_i64(1),
                vec![1, 0],
            )],
        )
        .expect("polynomial should exist");

        let point = [F17::from_i64(2)];
        let error = polynomial.evaluate(&point).expect_err("arity should fail");

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
                q_term(1, 2, &[1, 0]),
                q_term(2, 3, &[0, 1]),
                q_term(1, 4, &[0, 0]),
            ],
        )
        .expect("polynomial should exist");

        let point = [q(3, 2), q(6, 5)];
        let value = polynomial.evaluate(&point).expect("evaluation should work");

        assert!(Q::eq(&value, &q(9, 5)));
    }

    #[test]
    fn monomial_reports_arity_and_total_degree() {
        let monomial = Monomial::new(vec![2, 0, 1]);

        assert_eq!(monomial.arity(), 3);
        assert_eq!(monomial.total_degree(), 3);
    }

    #[test]
    fn monomial_multiplication_adds_exponents_component_wise() {
        let lhs = Monomial::new(vec![1, 2, 0]);
        let rhs = Monomial::new(vec![3, 0, 4]);
        let product = lhs.mul(&rhs).expect("arities should match");

        assert_eq!(product.exponents(), vec![4, 2, 4]);
    }

    #[test]
    fn monomial_multiplication_rejects_different_arities() {
        let lhs = Monomial::new(vec![1, 2]);
        let rhs = Monomial::new(vec![3, 4, 5]);

        assert!(lhs.mul(&rhs).is_none());
    }

    #[test]
    fn multivariate_polynomial_rejects_incompatible_arities() {
        let error = MultivariatePolynomial::<F17>::new(
            2,
            vec![f17_term(3, &[1, 0]), f17_term(4, &[0, 1, 2])],
        )
        .expect_err("mixed arities should fail");

        assert_eq!(
            error,
            PolynomialError::MonomialArityMismatch {
                expected: 2,
                actual: 3,
            }
        );
    }

    #[test]
    fn multivariate_polynomial_normalizes_zero_terms_and_collisions() {
        let polynomial = MultivariatePolynomial::<F17>::new(
            2,
            vec![
                f17_term(3, &[1, 0]),
                f17_term(14, &[1, 0]),
                f17_term(0, &[0, 2]),
                f17_term(5, &[0, 0]),
            ],
        )
        .expect("polynomial should exist");

        let terms = polynomial.terms();
        assert_eq!(polynomial.arity(), 2);
        assert_eq!(terms.len(), 1);
        assert_f17_term(&terms[0], &[0, 0], 5);
    }

    #[test]
    fn multivariate_polynomial_can_represent_zero_with_explicit_arity() {
        let polynomial =
            MultivariatePolynomial::<F17>::new(3, Vec::new()).expect("zero polynomial is valid");

        assert_eq!(polynomial.arity(), 3);
        assert!(polynomial.is_empty());
        assert_eq!(polynomial.len(), 0);
        assert_eq!(polynomial.degree(), None);
        assert!(polynomial.leading_term().is_none());
    }

    #[test]
    fn multivariate_polynomial_addition_combines_like_monomials() {
        let lhs =
            MultivariatePolynomial::<F17>::new(2, vec![f17_term(3, &[1, 0]), f17_term(5, &[0, 1])])
                .expect("lhs should exist");
        let rhs = MultivariatePolynomial::<F17>::new(
            2,
            vec![f17_term(14, &[1, 0]), f17_term(1, &[2, 0])],
        )
        .expect("rhs should exist");
        let sum = lhs.add(&rhs).expect("arities match");

        let terms = sum.terms();
        assert_eq!(terms.len(), 2);
        assert_f17_term(&terms[0], &[0, 1], 5);
        assert_f17_term(&terms[1], &[2, 0], 1);
    }

    #[test]
    fn multivariate_polynomial_addition_rejects_different_arities() {
        let lhs = MultivariatePolynomial::<F17>::new(2, vec![f17_term(1, &[1, 0])])
            .expect("lhs should exist");
        let rhs = MultivariatePolynomial::<F17>::new(3, vec![f17_term(1, &[1, 0, 0])])
            .expect("rhs should exist");

        let error = lhs.add(&rhs).expect_err("different arities should fail");
        assert_eq!(
            error,
            PolynomialError::IncompatibleMultivariateArity {
                lhs: 2,
                rhs: 3,
                operation: "addition",
            }
        );
    }

    #[test]
    fn multivariate_polynomial_multiplication_uses_naive_term_products() {
        let lhs =
            MultivariatePolynomial::<F17>::new(2, vec![f17_term(2, &[1, 0]), f17_term(3, &[0, 1])])
                .expect("lhs should exist");
        let rhs =
            MultivariatePolynomial::<F17>::new(2, vec![f17_term(4, &[0, 0]), f17_term(5, &[1, 0])])
                .expect("rhs should exist");
        let product = lhs.mul(&rhs).expect("arities match");

        let terms = product.terms();
        assert_eq!(terms.len(), 4);
        assert_f17_term(&terms[0], &[0, 1], 12);
        assert_f17_term(&terms[1], &[1, 0], 8);
        assert_f17_term(&terms[2], &[1, 1], 15);
        assert_f17_term(&terms[3], &[2, 0], 10);
        assert_eq!(product.degree(), Some(2));
    }

    #[test]
    fn multivariate_polynomial_multiplication_rejects_different_arities() {
        let lhs = MultivariatePolynomial::<F17>::new(2, vec![f17_term(1, &[1, 0])])
            .expect("lhs should exist");
        let rhs = MultivariatePolynomial::<F17>::new(1, vec![f17_term(1, &[2])])
            .expect("rhs should exist");

        let error = lhs.mul(&rhs).expect_err("different arities should fail");
        assert_eq!(
            error,
            PolynomialError::IncompatibleMultivariateArity {
                lhs: 2,
                rhs: 1,
                operation: "multiplication",
            }
        );
    }

    #[test]
    fn multivariate_polynomial_addition_works_over_q_too() {
        let lhs =
            MultivariatePolynomial::<Q>::new(2, vec![q_term(1, 2, &[1, 0]), q_term(2, 3, &[0, 1])])
                .expect("lhs should exist");
        let rhs = MultivariatePolynomial::<Q>::new(
            2,
            vec![
                q_term(1, 3, &[1, 0]),
                q_term(-2, 3, &[0, 1]),
                q_term(5, 4, &[0, 0]),
            ],
        )
        .expect("rhs should exist");
        let sum = lhs.add(&rhs).expect("arities match");

        let terms = sum.terms();
        assert_eq!(terms.len(), 2);
        assert_q_term(&terms[0], &[0, 0], 5, 4);
        assert_q_term(&terms[1], &[1, 0], 5, 6);
    }

    #[test]
    fn multivariate_polynomial_multiplication_works_over_q_too() {
        let lhs =
            MultivariatePolynomial::<Q>::new(2, vec![q_term(1, 2, &[1, 0]), q_term(1, 3, &[0, 1])])
                .expect("lhs should exist");
        let rhs =
            MultivariatePolynomial::<Q>::new(2, vec![q_term(2, 5, &[0, 0]), q_term(3, 7, &[1, 0])])
                .expect("rhs should exist");
        let product = lhs.mul(&rhs).expect("arities match");

        let terms = product.terms();
        assert_eq!(terms.len(), 4);
        assert_q_term(&terms[0], &[0, 1], 2, 15);
        assert_q_term(&terms[1], &[1, 0], 1, 5);
        assert_q_term(&terms[2], &[1, 1], 1, 7);
        assert_q_term(&terms[3], &[2, 0], 3, 14);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]

        #[test]
        fn property_multivariate_polynomials_stay_normalized(
            polynomial in arb_multivariate_polynomial::<F17>(PolynomialStrategyConfig {
                arity: 2,
                max_terms: 5,
                max_exponent: 3,
                ..PolynomialStrategyConfig::default()
            }),
        ) {
            let terms = polynomial.terms();
            prop_assert_eq!(polynomial.arity(), 2);
            prop_assert!(terms.iter().all(|term| term.monomial().arity() == polynomial.arity()));
            prop_assert!(terms.iter().all(|term| !F17::is_zero(term.coefficient())));
            prop_assert!(terms.windows(2).all(|window| window[0].monomial() < window[1].monomial()));
        }

        #[test]
        fn property_multivariate_evaluation_respects_multiplication(
            left in arb_multivariate_polynomial::<F17>(PolynomialStrategyConfig {
                arity: 2,
                max_terms: 4,
                max_exponent: 3,
                ..PolynomialStrategyConfig::default()
            }),
            right in arb_multivariate_polynomial::<F17>(PolynomialStrategyConfig {
                arity: 2,
                max_terms: 4,
                max_exponent: 3,
                ..PolynomialStrategyConfig::default()
            }),
            x0 in arb_fp_elem::<crate::fields::Fp17>(),
            x1 in arb_fp_elem::<crate::fields::Fp17>(),
        ) {
            let point = [x0, x1];
            let product = left.mul(&right).expect("matching arities should multiply");
            let product_value = product.evaluate(&point).expect("product should evaluate");
            let left_value = left.evaluate(&point).expect("left should evaluate");
            let right_value = right.evaluate(&point).expect("right should evaluate");

            prop_assert!(F17::eq(&product_value, &F17::mul(&left_value, &right_value)));
        }
    }
}
