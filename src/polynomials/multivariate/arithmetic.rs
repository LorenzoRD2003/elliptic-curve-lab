use crate::fields::traits::*;
use crate::polynomials::{MultivariatePolynomial, PolynomialError, multivariate::MultivariateTerm};

impl<F: Field> MultivariatePolynomial<F> {
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
