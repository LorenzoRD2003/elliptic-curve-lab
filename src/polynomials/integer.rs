use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

use crate::numerics::pow_bigint_usize;

#[derive(Clone, Debug, PartialEq, Eq)]
struct IntegerPolynomialTerm {
    degree: usize,
    coefficient: BigInt,
}

/// Sparse univariate polynomial with integer coefficients.
///
/// Terms are stored in strictly increasing degree order, and zero coefficients
/// are omitted. For example, `x² − 25` is stored as the two terms
/// `(0, −25)` and `(2, 1)`.
///
/// The constructor still accepts dense coefficients in ascending-degree order
/// because that is convenient for tests and small examples, but the internal
/// representation is sparse.
///
/// This is intentionally not wired into the field-oriented
/// [`crate::polynomials::DensePolynomial`] hierarchy. Integer polynomials live
/// in `ℤ[x]`, and `ℤ` is not a field, so operations that require division by an
/// arbitrary leading coefficient would be dishonest here.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IntegerPolynomial {
    terms: Vec<IntegerPolynomialTerm>,
}

impl IntegerPolynomial {
    /// Builds an integer polynomial from dense coefficients.
    ///
    /// Zero coefficients are skipped, so `[5, 0, 0]` stores one term and
    /// `[0, 0]` stores none.
    ///
    /// Complexity: `Θ(n)` in the number of supplied coefficients.
    pub(crate) fn new(coefficients: Vec<BigInt>) -> Self {
        Self {
            terms: coefficients
                .into_iter()
                .enumerate()
                .filter_map(|(degree, coefficient)| {
                    if coefficient.is_zero() {
                        None
                    } else {
                        Some(IntegerPolynomialTerm {
                            degree,
                            coefficient,
                        })
                    }
                })
                .collect(),
        }
    }

    /// Builds the candidate polynomial `x^q − N`.
    ///
    /// Complexity: `Θ(1)` terms are allocated for `q > 0`; the `q = 0` case is
    /// normalized as the constant polynomial `1 − N`.
    pub(crate) fn x_power_minus_constant(constant: &BigUint, exponent: u32) -> Option<Self> {
        let exponent = usize::try_from(exponent).ok()?;
        let constant = BigInt::from(constant.clone());
        if exponent == 0 {
            return Some(Self::new(vec![BigInt::one() - constant]));
        }

        Some(Self {
            terms: vec![
                IntegerPolynomialTerm {
                    degree: 0,
                    coefficient: -constant,
                },
                IntegerPolynomialTerm {
                    degree: exponent,
                    coefficient: BigInt::one(),
                },
            ],
        })
    }

    /// Materializes dense coefficients in ascending-degree order.
    ///
    /// This is mainly a compatibility bridge for lower-level Hensel traces that
    /// still record dense coefficient vectors.
    ///
    /// Complexity: `Θ(d + s)`, where `d = deg f` and `s` is the number of stored
    /// non-zero terms.
    pub(crate) fn to_dense_coefficients(&self) -> Vec<BigInt> {
        let Some(last_term) = self.terms.last() else {
            return Vec::new();
        };

        let mut coefficients = vec![BigInt::zero(); last_term.degree + 1];
        for term in &self.terms {
            coefficients[term.degree] = term.coefficient.clone();
        }
        coefficients
    }

    /// Returns whether this is the zero polynomial.
    pub(crate) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the degree of the polynomial when it is non-zero.
    pub(crate) fn degree(&self) -> Option<usize> {
        self.terms.last().map(|term| term.degree)
    }

    /// Returns whether this polynomial is constant.
    pub(crate) fn is_constant(&self) -> bool {
        self.degree() == Some(0)
    }

    /// Evaluates the polynomial at an integer point using sparse powers.
    ///
    /// Complexity: `Θ(s·log d)` exact integer multiplications for `s` stored
    /// terms and degree `d`, plus `Θ(s)` additions.
    pub(crate) fn evaluate(&self, x: &BigInt) -> BigInt {
        self.terms
            .iter()
            .map(|term| &term.coefficient * pow_bigint_usize(x, term.degree))
            .sum()
    }

    /// Evaluates the formal derivative at an integer point.
    ///
    /// Complexity: `Θ(s·log d)` exact integer multiplications for `s` stored
    /// terms and degree `d`, plus `Θ(s)` additions.
    pub(crate) fn evaluate_derivative(&self, x: &BigInt) -> BigInt {
        self.terms
            .iter()
            .filter(|term| term.degree > 0)
            .map(|term| {
                &term.coefficient * BigInt::from(term.degree) * pow_bigint_usize(x, term.degree - 1)
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::IntegerPolynomial;

    fn bi(value: i64) -> BigInt {
        BigInt::from(value)
    }

    #[test]
    fn dense_constructor_stores_sparse_terms_but_recovers_dense_coefficients() {
        let polynomial = IntegerPolynomial::new(vec![bi(-25), bi(0), bi(1)]);

        assert_eq!(polynomial.degree(), Some(2));
        assert!(!polynomial.is_constant());
        assert_eq!(
            polynomial.to_dense_coefficients(),
            vec![bi(-25), bi(0), bi(1)]
        );
    }

    #[test]
    fn sparse_evaluation_ignores_zero_gaps() {
        let polynomial = IntegerPolynomial::new(vec![bi(-25), bi(0), bi(1)]);

        assert_eq!(polynomial.evaluate(&bi(5)), bi(0));
        assert_eq!(polynomial.evaluate_derivative(&bi(5)), bi(10));
    }
}
