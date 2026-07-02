#[cfg(test)]
#[allow(clippy::module_inception)]
mod tests {
    use crate::fields::traits::*;

    use proptest::prelude::*;

    use crate::fields::{Q, traits::PthRootExtraction};
    use crate::polynomials::{
        DensePolynomial, SparsePolynomial, sparse::SparsePolynomialTerm,
        traits::UnivariatePolynomial,
    };
    use crate::proptest_support::{
        config::PolynomialStrategyConfig, fields::arb_fp_elem, polynomials::arb_sparse_polynomial,
    };

    type F17 = crate::fields::Fp17;

    crate::fields::extension_field::define_fp_quadratic_extension!(
        spec: F17Sqrt3SparsePthRootSpec,
        field: F17Sqrt3SparsePthRoot,
        base: F17,
        non_residue: 3,
        name: "F17(sqrt(3)) for sparse polynomial p-th-root tests",
    );

    fn f17_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::from_i64(coefficient),
            degree,
        }
    }

    fn q_term(numerator: i64, denominator: i64, degree: usize) -> SparsePolynomialTerm<Q> {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        SparsePolynomialTerm {
            coefficient: Q::div(&numerator, &denominator).expect("denominator should be non-zero"),
            degree,
        }
    }

    fn q(numerator: i64, denominator: i64) -> <Q as crate::fields::traits::Field>::Elem {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        Q::div(&numerator, &denominator).expect("denominator should be non-zero")
    }

    #[test]
    fn evaluate_sparse_matches_dense_over_f17() {
        let dense = DensePolynomial::<F17>::new(vec![
            F17::from_i64(3),
            F17::zero(),
            F17::from_i64(5),
            F17::from_i64(1),
        ]);
        let sparse =
            SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(5, 2), f17_term(1, 3)]);
        let point = F17::from_i64(2);

        let dense_value = dense
            .evaluate(&point)
            .expect("dense evaluation should work");
        let sparse_value = sparse
            .evaluate(&point)
            .expect("sparse evaluation should work");

        assert!(F17::eq(&dense_value, &sparse_value));
        assert!(F17::eq(&dense_value, &F17::from_i64(14)));
    }

    #[test]
    fn evaluate_sparse_works_over_q() {
        let polynomial = SparsePolynomial::<Q>::new(vec![q_term(1, 2, 0), q_term(2, 3, 2)]);
        let value = polynomial
            .evaluate(&q(3, 2))
            .expect("evaluation should work");

        assert!(Q::eq(&value, &q(2, 1)));
    }

    #[test]
    fn sparse_polynomial_normalizes_zero_terms_and_degree_collisions() {
        let polynomial = SparsePolynomial::<F17>::new(vec![
            f17_term(3, 2),
            f17_term(0, 5),
            f17_term(14, 2),
            f17_term(4, 0),
        ]);

        let terms = polynomial.terms();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].degree, 0);
        assert!(F17::eq(&terms[0].coefficient, &F17::from_i64(4)));
        assert!(polynomial.leading_term().is_some());
    }

    #[test]
    fn sparse_polynomial_cancels_terms_that_sum_to_zero() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(3, 1), f17_term(14, 1)]);

        assert!(polynomial.is_empty());
        assert_eq!(polynomial.degree(), None);
        assert_eq!(polynomial.leading_term(), None);
    }

    #[test]
    fn sparse_polynomial_preserves_sorted_non_zero_terms() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(5, 3), f17_term(2, 0)]);

        let terms = polynomial.terms();
        assert_eq!(terms.len(), 2);
        assert_eq!(terms[0].degree, 0);
        assert_eq!(terms[1].degree, 3);
        assert_eq!(polynomial.degree(), Some(3));
        assert_eq!(polynomial.leading_term().map(|term| term.degree), Some(3));
        assert!(F17::eq(
            polynomial
                .leading_coefficient()
                .expect("leading coefficient"),
            &F17::from_i64(5)
        ));
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::from_i64(2)
        ));
    }

    #[test]
    fn sparse_polynomial_addition_combines_matching_degrees() {
        let lhs = SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(5, 2)]);
        let rhs = SparsePolynomial::<F17>::new(vec![f17_term(14, 0), f17_term(1, 1)]);
        let sum = lhs.add(&rhs);

        let terms = sum.terms();
        assert_eq!(terms.len(), 2);
        assert_eq!(terms[0].degree, 1);
        assert!(F17::eq(&terms[0].coefficient, &F17::from_i64(1)));
        assert_eq!(terms[1].degree, 2);
        assert!(F17::eq(&terms[1].coefficient, &F17::from_i64(5)));
    }

    #[test]
    fn sparse_polynomial_constant_constructor_is_canonical() {
        let polynomial = SparsePolynomial::<F17>::constant(F17::from_i64(9));

        assert_eq!(polynomial.terms().len(), 1);
        assert_eq!(polynomial.terms()[0].degree, 0);
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::from_i64(9)
        ));

        let zero = SparsePolynomial::<F17>::constant(F17::zero());
        assert!(zero.is_zero());
        assert!(zero.is_empty());
    }

    #[test]
    fn sparse_polynomial_negation_and_subtraction_work_over_f17() {
        let lhs = SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(5, 2)]);
        let rhs = SparsePolynomial::<F17>::new(vec![f17_term(15, 0), f17_term(1, 1)]);

        let neg_rhs = rhs.neg();
        let expected_neg_rhs = SparsePolynomial::<F17>::new(vec![f17_term(2, 0), f17_term(16, 1)]);
        assert_eq!(neg_rhs, expected_neg_rhs);

        let difference = lhs.sub(&rhs);
        let expected_difference =
            SparsePolynomial::<F17>::new(vec![f17_term(5, 0), f17_term(16, 1), f17_term(5, 2)]);
        assert_eq!(difference, expected_difference);
    }

    #[test]
    fn sparse_polynomial_scale_multiplies_every_stored_coefficient() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(5, 2)]);
        let scaled = polynomial.scale(&F17::from_i64(4));

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(12, 0), f17_term(3, 2)]);
        assert_eq!(scaled, expected);
    }

    #[test]
    fn sparse_polynomial_derivative_drops_constant_terms_and_lowers_degrees() {
        let polynomial =
            SparsePolynomial::<F17>::new(vec![f17_term(4, 0), f17_term(3, 1), f17_term(5, 3)]);

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(15, 2)]);
        assert_eq!(polynomial.derivative(), expected);
    }

    #[test]
    fn sparse_polynomial_derivative_can_cancel_in_positive_characteristic() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(1, 17)]);

        assert!(polynomial.derivative().is_zero());
        assert_eq!(polynomial.derivative().terms(), &[]);
    }

    #[test]
    fn sparse_polynomial_derivative_works_over_q_too() {
        let polynomial =
            SparsePolynomial::<Q>::new(vec![q_term(1, 2, 0), q_term(2, 3, 1), q_term(3, 4, 2)]);

        let expected = SparsePolynomial::<Q>::new(vec![q_term(2, 3, 0), q_term(3, 2, 1)]);
        assert_eq!(polynomial.derivative(), expected);
    }

    #[test]
    fn sparse_polynomial_pth_root_over_prime_field_recovers_expected_terms() {
        let polynomial =
            SparsePolynomial::<F17>::new(vec![f17_term(4, 0), f17_term(9, 17), f17_term(3, 34)]);

        let root = polynomial
            .pth_root()
            .expect("all non-zero term degrees are divisible by the characteristic");

        let expected =
            SparsePolynomial::<F17>::new(vec![f17_term(4, 0), f17_term(9, 1), f17_term(3, 2)]);
        assert_eq!(root, expected);
        assert!(polynomial.has_pth_root());
    }

    #[test]
    fn sparse_polynomial_pth_root_rejects_non_multiple_degree_terms() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(1, 1), f17_term(3, 17)]);

        assert_eq!(polynomial.pth_root(), None);
        assert!(!polynomial.has_pth_root());
    }

    #[test]
    fn sparse_polynomial_pth_root_handles_zero_polynomial() {
        let polynomial = SparsePolynomial::<F17>::new(Vec::new());

        assert_eq!(
            polynomial.pth_root(),
            Some(SparsePolynomial::<F17>::new(Vec::new()))
        );
        assert!(polynomial.has_pth_root());
    }

    #[test]
    fn sparse_polynomial_pth_root_uses_extension_field_coefficient_roots() {
        let generator = F17Sqrt3SparsePthRoot::element(vec![F17::zero(), F17::one()]);
        let expected = SparsePolynomial::<F17Sqrt3SparsePthRoot>::new(vec![
            SparsePolynomialTerm {
                coefficient: generator.clone(),
                degree: 0,
            },
            SparsePolynomialTerm {
                coefficient: F17Sqrt3SparsePthRoot::one(),
                degree: 1,
            },
        ]);

        let polynomial = SparsePolynomial::<F17Sqrt3SparsePthRoot>::new(vec![
            SparsePolynomialTerm {
                coefficient: F17Sqrt3SparsePthRoot::pow(
                    &generator,
                    &F17Sqrt3SparsePthRoot::characteristic()
                        .to_positive_biguint()
                        .expect("finite fields should have positive characteristic"),
                ),
                degree: 0,
            },
            SparsePolynomialTerm {
                coefficient: F17Sqrt3SparsePthRoot::one(),
                degree: 17,
            },
        ]);

        assert!(polynomial.pth_root() == Some(expected));
        assert!(
            polynomial
                .constant_term()
                .expect("constant term should be present")
                .pth_root()
                == Some(generator)
        );
    }

    #[test]
    fn sparse_polynomial_gcd_returns_a_monic_common_divisor() {
        let lhs =
            SparsePolynomial::<F17>::new(vec![f17_term(2, 0), f17_term(3, 1), f17_term(1, 2)]);
        let rhs = SparsePolynomial::<F17>::new(vec![
            f17_term(1, 0),
            f17_term(3, 1),
            f17_term(3, 2),
            f17_term(1, 3),
        ]);

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(1, 0), f17_term(1, 1)]);
        let gcd = lhs.gcd(&rhs);
        assert_eq!(gcd, expected);
        assert!(gcd.is_monic());
    }

    #[test]
    fn sparse_polynomial_gcd_handles_zero_inputs() {
        let zero = SparsePolynomial::<F17>::new(Vec::new());
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(2, 0), f17_term(4, 1)]);

        assert!(zero.gcd(&zero).is_zero());
        assert_eq!(
            zero.gcd(&polynomial),
            SparsePolynomial::<F17>::new(vec![f17_term(9, 0), f17_term(1, 1)])
        );
        assert_eq!(
            polynomial.gcd(&zero),
            SparsePolynomial::<F17>::new(vec![f17_term(9, 0), f17_term(1, 1)])
        );
    }

    #[test]
    fn sparse_polynomial_multiplication_uses_naive_term_products() {
        let lhs = SparsePolynomial::<F17>::new(vec![f17_term(2, 0), f17_term(3, 1)]);
        let rhs = SparsePolynomial::<F17>::new(vec![f17_term(4, 0), f17_term(5, 2)]);
        let product = lhs.mul(&rhs);

        let terms = product.terms();
        assert_eq!(terms.len(), 4);
        assert_eq!(terms[0].degree, 0);
        assert!(F17::eq(&terms[0].coefficient, &F17::from_i64(8)));
        assert_eq!(terms[1].degree, 1);
        assert!(F17::eq(&terms[1].coefficient, &F17::from_i64(12)));
        assert_eq!(terms[2].degree, 2);
        assert!(F17::eq(&terms[2].coefficient, &F17::from_i64(10)));
        assert_eq!(terms[3].degree, 3);
        assert!(F17::eq(&terms[3].coefficient, &F17::from_i64(15)));
    }

    #[test]
    fn sparse_polynomial_addition_works_over_q_too() {
        let lhs = SparsePolynomial::<Q>::new(vec![q_term(1, 2, 0), q_term(2, 3, 2)]);
        let rhs = SparsePolynomial::<Q>::new(vec![q_term(1, 3, 0), q_term(-2, 3, 2)]);
        let sum = lhs.add(&rhs);

        let terms = sum.terms();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].degree, 0);
        let expected = Q::div(&Q::from_i64(5), &Q::from_i64(6)).unwrap();
        assert!(Q::eq(&terms[0].coefficient, &expected));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(36))]

        #[test]
        fn property_sparse_polynomials_stay_canonical(
            polynomial in arb_sparse_polynomial::<F17>(PolynomialStrategyConfig {
                max_terms: 6,
                max_degree: 6,
                ..PolynomialStrategyConfig::default()
            }),
        ) {
            let terms = polynomial.terms();
            prop_assert!(terms.iter().all(|term| !F17::is_zero(&term.coefficient)));
            prop_assert!(terms.windows(2).all(|window| window[0].degree < window[1].degree));
        }

        #[test]
        fn property_sparse_additive_inverse_cancels(
            polynomial in arb_sparse_polynomial::<F17>(PolynomialStrategyConfig {
                max_terms: 6,
                max_degree: 6,
                ..PolynomialStrategyConfig::default()
            }),
            scalar in arb_fp_elem::<crate::fields::Fp17>(),
        ) {
            let scaled = polynomial.scale(&scalar);
            prop_assert!(scaled.add(&scaled.neg()).is_zero());
            prop_assert_eq!(polynomial.sub(&polynomial), SparsePolynomial::<F17>::new(Vec::new()));
        }
    }

    #[test]
    fn sparse_polynomial_multiplication_works_over_q_too() {
        let lhs = SparsePolynomial::<Q>::new(vec![q_term(1, 2, 0), q_term(1, 3, 1)]);
        let rhs = SparsePolynomial::<Q>::new(vec![q_term(2, 5, 0), q_term(3, 7, 1)]);
        let product = lhs.mul(&rhs);

        let terms = product.terms();
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].degree, 0);
        assert!(Q::eq(
            &terms[0].coefficient,
            &Q::div(&Q::from_i64(1), &Q::from_i64(5)).unwrap()
        ));
        assert_eq!(terms[1].degree, 1);
        assert!(Q::eq(
            &terms[1].coefficient,
            &Q::div(&Q::from_i64(73), &Q::from_i64(210)).unwrap()
        ));
        assert_eq!(terms[2].degree, 2);
        assert!(Q::eq(
            &terms[2].coefficient,
            &Q::div(&Q::from_i64(1), &Q::from_i64(7)).unwrap()
        ));
    }

    #[test]
    fn sparse_polynomial_multiplication_preserves_empty_zero_representation() {
        let lhs = SparsePolynomial::<F17>::new(Vec::new());
        let rhs = SparsePolynomial::<F17>::new(vec![f17_term(1, 0), f17_term(2, 3)]);
        let product = lhs.mul(&rhs);

        assert!(product.is_zero());
        assert!(product.is_empty());
        assert_eq!(product.terms(), &[]);
    }

    #[test]
    fn dense_to_sparse_conversion_discards_zero_coefficients() {
        let dense = DensePolynomial::<F17>::new(vec![
            F17::from_i64(3),
            F17::zero(),
            F17::from_i64(5),
            F17::zero(),
            F17::from_i64(1),
        ]);
        let sparse = SparsePolynomial::<F17>::from(dense);

        let terms = sparse.terms();
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].degree, 0);
        assert!(F17::eq(&terms[0].coefficient, &F17::from_i64(3)));
        assert_eq!(terms[1].degree, 2);
        assert!(F17::eq(&terms[1].coefficient, &F17::from_i64(5)));
        assert_eq!(terms[2].degree, 4);
        assert!(F17::eq(&terms[2].coefficient, &F17::from_i64(1)));
    }

    #[test]
    fn dense_to_sparse_conversion_preserves_zero_polynomial() {
        let dense = DensePolynomial::<F17>::new(Vec::new());
        let sparse = SparsePolynomial::<F17>::from(dense);

        assert!(sparse.is_zero());
        assert!(sparse.is_empty());
        assert_eq!(sparse.terms(), &[]);
    }

    fn generic_scale<P>(polynomial: &P, scalar: &<F17 as crate::fields::traits::Field>::Elem) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        polynomial.scale(scalar)
    }

    fn generic_derivative<P>(polynomial: &P) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        polynomial.derivative()
    }

    fn generic_gcd<P>(lhs: &P, rhs: &P) -> P
    where
        P: UnivariatePolynomial<F17>,
    {
        lhs.gcd(rhs)
    }

    #[test]
    fn sparse_polynomial_implements_univariate_trait() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(3, 0), f17_term(5, 2)]);
        let scaled = generic_scale(&polynomial, &F17::from_i64(4));

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(12, 0), f17_term(3, 2)]);
        assert_eq!(scaled, expected);
        assert!(SparsePolynomial::<F17>::constant(F17::one()).is_monic());
    }

    #[test]
    fn sparse_polynomial_trait_derivative_uses_shared_surface() {
        let polynomial = SparsePolynomial::<F17>::new(vec![f17_term(6, 1), f17_term(5, 2)]);

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(6, 0), f17_term(10, 1)]);
        assert_eq!(generic_derivative(&polynomial), expected);
    }

    #[test]
    fn sparse_polynomial_trait_gcd_uses_shared_surface() {
        let lhs =
            SparsePolynomial::<F17>::new(vec![f17_term(2, 0), f17_term(3, 1), f17_term(1, 2)]);
        let rhs = SparsePolynomial::<F17>::new(vec![
            f17_term(1, 0),
            f17_term(3, 1),
            f17_term(3, 2),
            f17_term(1, 3),
        ]);

        let expected = SparsePolynomial::<F17>::new(vec![f17_term(1, 0), f17_term(1, 1)]);
        assert_eq!(generic_gcd(&lhs, &rhs), expected);
    }
}
