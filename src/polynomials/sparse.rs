use std::collections::BTreeMap;

use crate::fields::Field;
use crate::polynomials::DensePolynomial;
use crate::polynomials::UnivariatePolynomial;

/// One non-zero term of a sparse polynomial over a field `F`.
///
/// The term represents
///
/// `coefficient * x^degree`
///
/// As in the rest of the current `polynomials` module, the project is
/// intentionally specializing polynomial arithmetic to coefficient fields for
/// now, even though the mathematically more general story would allow arbitrary
/// coefficient rings.
#[derive(Debug)]
pub struct SparsePolynomialTerm<F: Field> {
    /// Non-zero coefficient of the term.
    pub coefficient: F::Elem,
    /// Degree of the monomial attached to the coefficient.
    pub degree: usize,
}

/// Sparse polynomial over a field `F`.
///
/// This representation stores only explicit non-zero terms and keeps them
/// normalized in ascending degree order.
///
/// The constructor performs a light but useful normalization pass:
///
/// - zero coefficients are discarded
/// - repeated degrees are combined by addition
/// - terms are stored in ascending degree order
///
/// This is intentionally stricter than the current dense representation, where
/// trailing-zero canonicalization is still deferred. For sparse polynomials,
/// normalization is part of what makes the representation pedagogically useful.
#[derive(Debug)]
pub struct SparsePolynomial<F: Field> {
    terms: Vec<SparsePolynomialTerm<F>>,
}

impl<F: Field> Clone for SparsePolynomialTerm<F> {
    fn clone(&self) -> Self {
        Self {
            coefficient: self.coefficient.clone(),
            degree: self.degree,
        }
    }
}

impl<F: Field> PartialEq for SparsePolynomialTerm<F> {
    fn eq(&self, other: &Self) -> bool {
        self.degree == other.degree && F::eq(&self.coefficient, &other.coefficient)
    }
}

impl<F: Field> Clone for SparsePolynomial<F> {
    fn clone(&self) -> Self {
        Self {
            terms: self.terms.clone(),
        }
    }
}

impl<F: Field> PartialEq for SparsePolynomial<F> {
    fn eq(&self, other: &Self) -> bool {
        self.terms == other.terms
    }
}

impl<F: Field> SparsePolynomial<F> {
    /// Builds a sparse polynomial from explicit terms.
    ///
    /// Terms with the same degree are added together. Terms whose coefficient
    /// is zero are removed from the final representation.
    pub fn new(terms: Vec<SparsePolynomialTerm<F>>) -> Self {
        Self {
            terms: Self::normalize_terms(terms),
        }
    }

    /// Builds the constant sparse polynomial with the given value.
    ///
    /// If the constant is zero, the result is the canonical zero polynomial
    /// with no stored terms.
    pub fn constant(value: F::Elem) -> Self {
        Self::new(vec![SparsePolynomialTerm {
            coefficient: value,
            degree: 0,
        }])
    }

    /// Normalizes a sparse term list into canonical storage order.
    ///
    /// The normalization pass:
    ///
    /// - discards zero coefficients
    /// - combines repeated degrees by addition
    /// - stores the surviving terms in ascending degree order
    fn normalize_terms(terms: Vec<SparsePolynomialTerm<F>>) -> Vec<SparsePolynomialTerm<F>> {
        let mut by_degree: BTreeMap<usize, F::Elem> = BTreeMap::new();

        for term in terms {
            let updated = if let Some(existing) = by_degree.remove(&term.degree) {
                F::add(&existing, &term.coefficient)
            } else {
                term.coefficient
            };

            if !F::is_zero(&updated) {
                by_degree.insert(term.degree, updated);
            }
        }

        by_degree
            .into_iter()
            .map(|(degree, coefficient)| SparsePolynomialTerm {
                coefficient,
                degree,
            })
            .collect()
    }

    /// Returns the stored terms in ascending degree order.
    pub fn terms(&self) -> &[SparsePolynomialTerm<F>] {
        &self.terms
    }

    /// Returns the number of stored non-zero terms.
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns the maximum stored degree.
    pub fn degree(&self) -> Option<usize> {
        self.terms.last().map(|term| term.degree)
    }

    /// Returns the leading term, if any.
    pub fn leading_term(&self) -> Option<&SparsePolynomialTerm<F>> {
        self.terms.last()
    }

    /// Returns the leading coefficient, if any.
    pub fn leading_coefficient(&self) -> Option<&F::Elem> {
        self.leading_term().map(|term| &term.coefficient)
    }

    /// Returns the constant term, if any.
    ///
    /// This is the coefficient attached to degree `0`. If no constant term is
    /// stored, the sparse representation treats that as an implicit zero and
    /// this method returns `None`.
    pub fn constant_term(&self) -> Option<&F::Elem> {
        self.terms
            .first()
            .filter(|term| term.degree == 0)
            .map(|term| &term.coefficient)
    }

    /// Returns whether the polynomial is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns whether the sparse representation stores no non-zero terms.
    ///
    /// This remains a useful storage-oriented query for sparse polynomials even
    /// though the algebraic query [`SparsePolynomial::is_zero`] also exists.
    pub fn is_empty(&self) -> bool {
        self.is_zero()
    }

    /// Adds two sparse polynomials.
    ///
    /// The implementation concatenates the term lists and lets the constructor
    /// perform normalization via [`SparsePolynomial::normalize_terms`]. This
    /// is simple and sufficient for the current educational stage of the
    /// project.
    pub fn add(&self, rhs: &Self) -> Self {
        let mut terms = Vec::with_capacity(self.len() + rhs.len());
        terms.extend(self.terms.iter().cloned());
        terms.extend(rhs.terms.iter().cloned());
        Self::new(terms)
    }

    /// Negates every stored coefficient of the polynomial.
    pub fn neg(&self) -> Self {
        Self::new(
            self.terms
                .iter()
                .map(|term| SparsePolynomialTerm {
                    coefficient: F::neg(&term.coefficient),
                    degree: term.degree,
                })
                .collect(),
        )
    }

    /// Subtracts two sparse polynomials.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// Multiplies every stored coefficient by the same field element.
    pub fn scale(&self, scalar: &F::Elem) -> Self {
        Self::new(
            self.terms
                .iter()
                .map(|term| SparsePolynomialTerm {
                    coefficient: F::mul(&term.coefficient, scalar),
                    degree: term.degree,
                })
                .collect(),
        )
    }

    /// Multiplies two sparse polynomials using the naive pairwise product.
    ///
    /// Every term on the left is multiplied by every term on the right. The
    /// constructor is then used to normalize the resulting term list via
    /// [`SparsePolynomial::normalize_terms`].
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_empty() || rhs.is_empty() {
            return Self::new(Vec::new());
        }

        let mut terms = Vec::with_capacity(self.len() * rhs.len());

        for lhs_term in &self.terms {
            for rhs_term in &rhs.terms {
                terms.push(SparsePolynomialTerm {
                    coefficient: F::mul(&lhs_term.coefficient, &rhs_term.coefficient),
                    degree: lhs_term.degree + rhs_term.degree,
                });
            }
        }

        Self::new(terms)
    }

    /// Returns the formal derivative of the polynomial.
    ///
    /// Each term `a*x^d` with `d > 0` becomes `d*a*x^(d-1)`. Constant terms
    /// disappear, and any coefficient that becomes zero in the base field is
    /// removed by the sparse constructor's normalization.
    pub fn derivative(&self) -> Self {
        let terms = self
            .terms
            .iter()
            .filter(|term| term.degree > 0)
            .map(|term| SparsePolynomialTerm {
                coefficient: F::mul(
                    &term.coefficient,
                    &F::elem_from_u64(
                        u64::try_from(term.degree)
                            .expect("sparse polynomial degree should fit in u64"),
                    ),
                ),
                degree: term.degree - 1,
            })
            .collect();

        Self::new(terms)
    }

    /// Computes a greatest common divisor by delegating to the current dense
    /// Euclidean algorithm and converting the monic result back to sparse
    /// storage.
    ///
    /// This keeps the shared algebraic answer available for sparse
    /// polynomials without yet duplicating polynomial division machinery in
    /// the sparse representation.
    pub fn gcd(&self, rhs: &Self) -> Self {
        let left_dense = DensePolynomial::<F>::from(self.clone());
        let right_dense = DensePolynomial::<F>::from(rhs.clone());
        let gcd_dense = left_dense.gcd(&right_dense);
        Self::from(gcd_dense)
    }
}

impl<F: Field> From<DensePolynomial<F>> for SparsePolynomial<F> {
    /// Converts a dense polynomial into its canonical sparse term list.
    ///
    /// Zero coefficients are discarded, so only the non-zero degrees survive
    /// in the sparse representation.
    fn from(polynomial: DensePolynomial<F>) -> Self {
        let terms = polynomial
            .coefficients()
            .iter()
            .enumerate()
            .filter_map(|(degree, coefficient)| {
                if F::is_zero(coefficient) {
                    None
                } else {
                    Some(SparsePolynomialTerm {
                        coefficient: coefficient.clone(),
                        degree,
                    })
                }
            })
            .collect();

        Self::new(terms)
    }
}

impl<F: Field> UnivariatePolynomial<F> for SparsePolynomial<F> {
    fn constant(value: F::Elem) -> Self {
        SparsePolynomial::constant(value)
    }

    fn degree(&self) -> Option<usize> {
        SparsePolynomial::degree(self)
    }

    fn leading_coefficient(&self) -> Option<&F::Elem> {
        SparsePolynomial::leading_coefficient(self)
    }

    fn constant_term(&self) -> Option<&F::Elem> {
        SparsePolynomial::constant_term(self)
    }

    fn is_zero(&self) -> bool {
        SparsePolynomial::is_zero(self)
    }

    fn add(&self, rhs: &Self) -> Self {
        SparsePolynomial::add(self, rhs)
    }

    fn neg(&self) -> Self {
        SparsePolynomial::neg(self)
    }

    fn sub(&self, rhs: &Self) -> Self {
        SparsePolynomial::sub(self, rhs)
    }

    fn scale(&self, scalar: &F::Elem) -> Self {
        SparsePolynomial::scale(self, scalar)
    }

    fn mul(&self, rhs: &Self) -> Self {
        SparsePolynomial::mul(self, rhs)
    }

    fn derivative(&self) -> Self {
        SparsePolynomial::derivative(self)
    }

    fn gcd(&self, rhs: &Self) -> Self {
        SparsePolynomial::gcd(self, rhs)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::fields::{Field, Fp, Q};
    use crate::polynomials::{DensePolynomial, UnivariatePolynomial};
    use crate::proptest_support::config::PolynomialStrategyConfig;
    use crate::proptest_support::fields::arb_fp_elem;
    use crate::proptest_support::polynomials::arb_sparse_polynomial;

    use crate::polynomials::{SparsePolynomial, SparsePolynomialTerm};

    type F17 = Fp<17>;

    fn f17_term(coefficient: u64, degree: usize) -> SparsePolynomialTerm<F17> {
        SparsePolynomialTerm {
            coefficient: F17::elem_from_u64(coefficient),
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
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(4)));
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
            &F17::elem_from_u64(5)
        ));
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::elem_from_u64(2)
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
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(1)));
        assert_eq!(terms[1].degree, 2);
        assert!(F17::eq(&terms[1].coefficient, &F17::elem_from_u64(5)));
    }

    #[test]
    fn sparse_polynomial_constant_constructor_is_canonical() {
        let polynomial = SparsePolynomial::<F17>::constant(F17::elem_from_u64(9));

        assert_eq!(polynomial.terms().len(), 1);
        assert_eq!(polynomial.terms()[0].degree, 0);
        assert!(F17::eq(
            polynomial.constant_term().expect("constant term"),
            &F17::elem_from_u64(9)
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
        let scaled = polynomial.scale(&F17::elem_from_u64(4));

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
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(8)));
        assert_eq!(terms[1].degree, 1);
        assert!(F17::eq(&terms[1].coefficient, &F17::elem_from_u64(12)));
        assert_eq!(terms[2].degree, 2);
        assert!(F17::eq(&terms[2].coefficient, &F17::elem_from_u64(10)));
        assert_eq!(terms[3].degree, 3);
        assert!(F17::eq(&terms[3].coefficient, &F17::elem_from_u64(15)));
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
            scalar in arb_fp_elem::<17>(),
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
            F17::elem_from_u64(3),
            F17::zero(),
            F17::elem_from_u64(5),
            F17::zero(),
            F17::elem_from_u64(1),
        ]);
        let sparse = SparsePolynomial::<F17>::from(dense);

        let terms = sparse.terms();
        assert_eq!(terms.len(), 3);
        assert_eq!(terms[0].degree, 0);
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(3)));
        assert_eq!(terms[1].degree, 2);
        assert!(F17::eq(&terms[1].coefficient, &F17::elem_from_u64(5)));
        assert_eq!(terms[2].degree, 4);
        assert!(F17::eq(&terms[2].coefficient, &F17::elem_from_u64(1)));
    }

    #[test]
    fn dense_to_sparse_conversion_preserves_zero_polynomial() {
        let dense = DensePolynomial::<F17>::new(Vec::new());
        let sparse = SparsePolynomial::<F17>::from(dense);

        assert!(sparse.is_zero());
        assert!(sparse.is_empty());
        assert_eq!(sparse.terms(), &[]);
    }

    fn generic_scale<P>(polynomial: &P, scalar: &<F17 as Field>::Elem) -> P
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
        let scaled = generic_scale(&polynomial, &F17::elem_from_u64(4));

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
