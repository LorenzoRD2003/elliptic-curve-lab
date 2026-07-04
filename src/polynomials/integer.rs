use num_bigint::{BigInt, BigUint};
use num_traits::{One, Signed, Zero};

use crate::numerics::{positive_divisors, pow_bigint_usize};

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
    #[cfg(test)]
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

    /// Returns a Cauchy bound for the absolute value of integer roots.
    ///
    /// For a nonconstant polynomial
    /// `f(x) = a_d x^d + a_{d-1}x^{d-1} + ⋯ + a_0`, every complex root
    /// satisfies
    ///
    /// `|x| ≤ 1 + max_{i < d} |a_i| / |a_d|`.
    ///
    /// This helper returns the integer ceiling of that bound.
    /// It returns `None` for constant polynomials.
    ///
    /// Complexity: `Θ(s)` integer operations for `s` stored non-zero terms.
    pub(crate) fn cauchy_integer_root_bound(&self) -> Option<BigUint> {
        let leading = self.terms.last()?;
        if leading.degree == 0 {
            return None;
        }

        let leading_abs = leading
            .coefficient
            .abs()
            .to_biguint()
            .expect("absolute leading coefficient should be nonnegative");
        let mut max_ratio_ceiling = BigUint::zero();

        for term in self
            .terms
            .iter()
            .filter(|term| term.degree < leading.degree)
        {
            let coefficient_abs = term
                .coefficient
                .abs()
                .to_biguint()
                .expect("absolute coefficient should be nonnegative");
            let ratio_ceiling = ceil_div_biguint(&coefficient_abs, &leading_abs);
            if ratio_ceiling > max_ratio_ceiling {
                max_ratio_ceiling = ratio_ceiling;
            }
        }

        Some(max_ratio_ceiling + BigUint::one())
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

    /// Finds integer roots using the rational-root theorem.
    ///
    /// This is intentionally a small exact helper for low-degree educational
    /// workflows. If the constant term vanishes, the root `0` is recorded and
    /// factors of `x` are stripped before enumerating divisors of the first
    /// non-zero constant term. The zero polynomial has infinitely many roots,
    /// so this finite-root helper returns an empty list for it.
    ///
    /// Complexity: `Θ(factor(c) + τ(c)·eval)`, where `c` is the first non-zero
    /// constant term after stripping powers of `x`, `τ(c)` is its number of
    /// positive divisors, and `eval` is the cost of one sparse evaluation.
    pub(crate) fn integer_roots_by_rational_root_test(&self) -> Vec<BigInt> {
        if self.is_zero() {
            return Vec::new();
        }

        let mut reduced = self.clone();
        let mut roots = Vec::new();

        while !reduced.is_zero() && reduced.constant_coefficient().is_zero() {
            if roots.is_empty() {
                roots.push(BigInt::zero());
            }
            reduced = reduced
                .divide_by_x()
                .expect("zero constant term should make division by x exact");
        }

        if reduced.is_zero() {
            roots.sort();
            return roots;
        }

        let constant = reduced.constant_coefficient().abs();
        let divisors = positive_divisors(
            &constant
                .to_biguint()
                .expect("absolute constant coefficient should be nonnegative"),
        );
        for divisor in divisors {
            let positive = BigInt::from(divisor);
            self.push_root_if_exact(&positive, &mut roots);
            if !positive.is_zero() {
                self.push_root_if_exact(&(-positive), &mut roots);
            }
        }

        roots.sort();
        roots.dedup();
        roots
    }

    fn constant_coefficient(&self) -> BigInt {
        self.terms
            .first()
            .filter(|term| term.degree == 0)
            .map(|term| term.coefficient.clone())
            .unwrap_or_else(BigInt::zero)
    }

    fn divide_by_x(&self) -> Option<Self> {
        if !self.constant_coefficient().is_zero() {
            return None;
        }

        Some(Self {
            terms: self
                .terms
                .iter()
                .filter(|term| term.degree > 0)
                .map(|term| IntegerPolynomialTerm {
                    degree: term.degree - 1,
                    coefficient: term.coefficient.clone(),
                })
                .collect(),
        })
    }

    fn push_root_if_exact(&self, candidate: &BigInt, roots: &mut Vec<BigInt>) {
        if self.evaluate(candidate).is_zero() {
            roots.push(candidate.clone());
        }
    }
}

fn ceil_div_biguint(numerator: &BigUint, denominator: &BigUint) -> BigUint {
    debug_assert!(!denominator.is_zero());
    if numerator.is_zero() {
        BigUint::zero()
    } else {
        (numerator + denominator - BigUint::one()) / denominator
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_bigint::BigUint;

    use super::IntegerPolynomial;

    fn bi(value: i64) -> BigInt {
        BigInt::from(value)
    }

    fn bu(value: u64) -> BigUint {
        BigUint::from(value)
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

    #[test]
    fn integer_roots_include_zero_after_stripping_x_factors() {
        let polynomial = IntegerPolynomial::new(vec![bi(0), bi(-1), bi(0), bi(1)]);

        assert_eq!(
            polynomial.integer_roots_by_rational_root_test(),
            vec![bi(-1), bi(0), bi(1)]
        );
    }

    #[test]
    fn integer_roots_use_signed_divisors_of_nonzero_constant() {
        let polynomial = IntegerPolynomial::new(vec![bi(-6), bi(11), bi(-6), bi(1)]);

        assert_eq!(
            polynomial.integer_roots_by_rational_root_test(),
            vec![bi(1), bi(2), bi(3)]
        );
    }

    #[test]
    fn cauchy_integer_root_bound_bounds_known_integer_roots() {
        let polynomial = IntegerPolynomial::new(vec![bi(-6), bi(11), bi(-6), bi(1)]);

        assert_eq!(polynomial.cauchy_integer_root_bound(), Some(bu(12)));
    }

    #[test]
    fn cauchy_integer_root_bound_handles_sparse_monic_polynomials() {
        let polynomial = IntegerPolynomial::new(vec![bi(-25), bi(0), bi(1)]);

        assert_eq!(polynomial.cauchy_integer_root_bound(), Some(bu(26)));
    }

    #[test]
    fn cauchy_integer_root_bound_rejects_zero_and_constant_polynomials() {
        assert_eq!(
            IntegerPolynomial::new(Vec::new()).cauchy_integer_root_bound(),
            None
        );
        assert_eq!(
            IntegerPolynomial::new(vec![bi(5)]).cauchy_integer_root_bound(),
            None
        );
    }
}
