use std::collections::BTreeMap;

use crate::fields::traits::Field;
use crate::polynomials::{PolynomialError, multivariate::Monomial};

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
    pub(super) arity: usize,
    pub(super) terms: Vec<MultivariateTerm<F>>,
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
    pub(super) fn normalize_terms(
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
}
