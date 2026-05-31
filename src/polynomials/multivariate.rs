use std::collections::BTreeMap;

use crate::fields::Field;
use crate::polynomials::PolynomialError;

/// Exponent vector of a multivariate monomial.
///
/// If the ambient variables are ordered as
///
/// `x_0, x_1, ..., x_{n-1}`
///
/// then the exponent vector `[e0, e1, ..., e_{n-1}]` represents
///
/// `x_0^e0 * x_1^e1 * ... * x_{n-1}^e_{n-1}`
///
/// For example, in arity `3`:
///
/// - `[2, 0, 1]` represents `x_0^2 * x_2`
/// - `[0, 1, 0]` represents `x_1`
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Monomial {
    /// Exponents of the ordered variable list.
    pub exponents: Vec<usize>,
}

impl Monomial {
    /// Builds a monomial from an explicit exponent vector.
    pub fn new(exponents: Vec<usize>) -> Self {
        Self { exponents }
    }

    /// Returns the ambient number of variables.
    pub fn arity(&self) -> usize {
        self.exponents.len()
    }

    /// Returns the total degree of the monomial.
    pub fn total_degree(&self) -> usize {
        self.exponents.iter().sum()
    }

    /// Multiplies two monomials of the same arity by adding exponents
    /// component-wise.
    pub fn mul(&self, rhs: &Self) -> Option<Self> {
        if self.arity() != rhs.arity() {
            return None;
        }

        let exponents = self
            .exponents
            .iter()
            .zip(&rhs.exponents)
            .map(|(lhs, rhs)| lhs + rhs)
            .collect();

        Some(Self { exponents })
    }
}

/// One non-zero term of a multivariate polynomial over a field `F`.
#[derive(Debug)]
pub struct MultivariateTerm<F: Field> {
    /// Non-zero coefficient of the term.
    pub coefficient: F::Elem,
    /// Monomial attached to the coefficient.
    pub monomial: Monomial,
}

impl<F: Field> Clone for MultivariateTerm<F> {
    fn clone(&self) -> Self {
        Self {
            coefficient: self.coefficient.clone(),
            monomial: self.monomial.clone(),
        }
    }
}

impl<F: Field> PartialEq for MultivariateTerm<F> {
    fn eq(&self, other: &Self) -> bool {
        self.monomial == other.monomial && F::eq(&self.coefficient, &other.coefficient)
    }
}

/// Multivariate polynomial over a field `F`.
///
/// This representation is intentionally specialized to coefficient fields in
/// the current stage of the project. The more general algebraic story would
/// allow coefficients from broader classes of rings, but the present goal is
/// to support polynomial arithmetic over the field types that already exist in
/// the repository.
///
/// The polynomial stores:
///
/// - an explicit ambient arity
/// - a normalized sparse list of non-zero terms
///
/// The constructor enforces these invariants:
///
/// - every stored monomial has the declared arity
/// - zero coefficients are discarded
/// - repeated monomials are combined by addition
/// - terms are stored in ascending monomial order
#[derive(Debug)]
pub struct MultivariatePolynomial<F: Field> {
    arity: usize,
    terms: Vec<MultivariateTerm<F>>,
}

impl<F: Field> Clone for MultivariatePolynomial<F> {
    fn clone(&self) -> Self {
        Self {
            arity: self.arity,
            terms: self.terms.clone(),
        }
    }
}

impl<F: Field> PartialEq for MultivariatePolynomial<F> {
    fn eq(&self, other: &Self) -> bool {
        self.arity == other.arity && self.terms == other.terms
    }
}

impl<F: Field> MultivariatePolynomial<F> {
    /// Builds a multivariate polynomial from explicit terms in a fixed ambient
    /// arity.
    pub fn new(arity: usize, terms: Vec<MultivariateTerm<F>>) -> Result<Self, PolynomialError> {
        let terms = Self::normalize_terms(arity, terms)?;

        Ok(Self { arity, terms })
    }

    /// Normalizes a multivariate term list for a fixed ambient arity.
    ///
    /// The normalization pass:
    ///
    /// - rejects monomials whose exponent vector does not match the declared
    ///   arity
    /// - discards zero coefficients
    /// - combines repeated monomials by addition
    /// - stores the surviving terms in ascending monomial order
    fn normalize_terms(
        arity: usize,
        terms: Vec<MultivariateTerm<F>>,
    ) -> Result<Vec<MultivariateTerm<F>>, PolynomialError> {
        let mut by_monomial: BTreeMap<Monomial, F::Elem> = BTreeMap::new();

        for term in terms {
            if term.monomial.arity() != arity {
                return Err(PolynomialError::MonomialArityMismatch {
                    expected: arity,
                    actual: term.monomial.arity(),
                });
            }

            let updated = if let Some(existing) = by_monomial.remove(&term.monomial) {
                F::add(&existing, &term.coefficient)
            } else {
                term.coefficient
            };

            if !F::is_zero(&updated) {
                by_monomial.insert(term.monomial, updated);
            }
        }

        Ok(by_monomial
            .into_iter()
            .map(|(monomial, coefficient)| MultivariateTerm {
                coefficient,
                monomial,
            })
            .collect())
    }

    /// Returns the ambient number of variables.
    pub fn arity(&self) -> usize {
        self.arity
    }

    /// Returns the stored non-zero terms in ascending monomial order.
    pub fn terms(&self) -> &[MultivariateTerm<F>] {
        &self.terms
    }

    /// Returns the number of stored non-zero terms.
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Returns whether the polynomial stores no non-zero terms.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the maximum total degree among the stored terms.
    pub fn degree(&self) -> Option<usize> {
        self.terms
            .iter()
            .map(|term| term.monomial.total_degree())
            .max()
    }

    /// Returns the last term in the stored monomial order, if any.
    pub fn leading_term(&self) -> Option<&MultivariateTerm<F>> {
        self.terms.last()
    }

    /// Adds two multivariate polynomials of the same arity.
    ///
    /// The concatenated term list is normalized through
    /// [`MultivariatePolynomial::normalize_terms`].
    pub fn add(&self, rhs: &Self) -> Result<Self, PolynomialError> {
        if self.arity != rhs.arity {
            return Err(PolynomialError::IncompatibleMultivariateArity {
                lhs: self.arity,
                rhs: rhs.arity,
                operation: "addition",
            });
        }

        let mut terms = Vec::with_capacity(self.len() + rhs.len());
        terms.extend(self.terms.iter().cloned());
        terms.extend(rhs.terms.iter().cloned());
        Self::new(self.arity, terms)
    }

    /// Multiplies two multivariate polynomials of the same arity.
    ///
    /// The naive pairwise product is normalized through
    /// [`MultivariatePolynomial::normalize_terms`].
    pub fn mul(&self, rhs: &Self) -> Result<Self, PolynomialError> {
        if self.arity != rhs.arity {
            return Err(PolynomialError::IncompatibleMultivariateArity {
                lhs: self.arity,
                rhs: rhs.arity,
                operation: "multiplication",
            });
        }

        if self.is_empty() || rhs.is_empty() {
            return Self::new(self.arity, Vec::new());
        }

        let mut terms = Vec::with_capacity(self.len() * rhs.len());

        for lhs_term in &self.terms {
            for rhs_term in &rhs.terms {
                let monomial = lhs_term
                    .monomial
                    .mul(&rhs_term.monomial)
                    .expect("polynomials with matching arity have compatible monomials");

                terms.push(MultivariateTerm {
                    coefficient: F::mul(&lhs_term.coefficient, &rhs_term.coefficient),
                    monomial,
                });
            }
        }

        Self::new(self.arity, terms)
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::{Field, Fp, Q};
    use crate::polynomials::PolynomialError;

    use super::{Monomial, MultivariatePolynomial, MultivariateTerm};

    type F17 = Fp<17>;

    fn f17_term(coefficient: u64, exponents: &[usize]) -> MultivariateTerm<F17> {
        MultivariateTerm {
            coefficient: F17::elem_from_u64(coefficient),
            monomial: Monomial::new(exponents.to_vec()),
        }
    }

    fn q_term(numerator: i64, denominator: i64, exponents: &[usize]) -> MultivariateTerm<Q> {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        MultivariateTerm {
            coefficient: Q::div(&numerator, &denominator).expect("denominator should be non-zero"),
            monomial: Monomial::new(exponents.to_vec()),
        }
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

        assert_eq!(product.exponents, vec![4, 2, 4]);
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
        assert_eq!(terms[0].monomial.exponents, vec![0, 0]);
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(5)));
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
        assert_eq!(terms[0].monomial.exponents, vec![0, 1]);
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(5)));
        assert_eq!(terms[1].monomial.exponents, vec![2, 0]);
        assert!(F17::eq(&terms[1].coefficient, &F17::elem_from_u64(1)));
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
        assert_eq!(terms[0].monomial.exponents, vec![0, 1]);
        assert!(F17::eq(&terms[0].coefficient, &F17::elem_from_u64(12)));
        assert_eq!(terms[1].monomial.exponents, vec![1, 0]);
        assert!(F17::eq(&terms[1].coefficient, &F17::elem_from_u64(8)));
        assert_eq!(terms[2].monomial.exponents, vec![1, 1]);
        assert!(F17::eq(&terms[2].coefficient, &F17::elem_from_u64(15)));
        assert_eq!(terms[3].monomial.exponents, vec![2, 0]);
        assert!(F17::eq(&terms[3].coefficient, &F17::elem_from_u64(10)));
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
        assert_eq!(terms[0].monomial.exponents, vec![0, 0]);
        assert!(Q::eq(
            &terms[0].coefficient,
            &Q::div(&Q::from_i64(5), &Q::from_i64(4)).unwrap()
        ));
        assert_eq!(terms[1].monomial.exponents, vec![1, 0]);
        assert!(Q::eq(
            &terms[1].coefficient,
            &Q::div(&Q::from_i64(5), &Q::from_i64(6)).unwrap()
        ));
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
        assert_eq!(terms[0].monomial.exponents, vec![0, 1]);
        assert!(Q::eq(
            &terms[0].coefficient,
            &Q::div(&Q::from_i64(2), &Q::from_i64(15)).unwrap()
        ));
        assert_eq!(terms[1].monomial.exponents, vec![1, 0]);
        assert!(Q::eq(
            &terms[1].coefficient,
            &Q::div(&Q::from_i64(1), &Q::from_i64(5)).unwrap()
        ));
        assert_eq!(terms[2].monomial.exponents, vec![1, 1]);
        assert!(Q::eq(
            &terms[2].coefficient,
            &Q::div(&Q::from_i64(1), &Q::from_i64(7)).unwrap()
        ));
        assert_eq!(terms[3].monomial.exponents, vec![2, 0]);
        assert!(Q::eq(
            &terms[3].coefficient,
            &Q::div(&Q::from_i64(3), &Q::from_i64(14)).unwrap()
        ));
    }
}
