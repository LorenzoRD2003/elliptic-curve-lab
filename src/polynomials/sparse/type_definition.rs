use crate::fields::traits::*;
use std::collections::BTreeMap;

/// One non-zero term of a sparse polynomial over a field `F`.
///
/// The term represents `coefficient * x^degree`
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
    pub(super) terms: Vec<SparsePolynomialTerm<F>>,
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
    pub(super) fn normalize_terms(
        terms: Vec<SparsePolynomialTerm<F>>,
    ) -> Vec<SparsePolynomialTerm<F>> {
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
}
