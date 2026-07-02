use crate::fields::traits::*;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

use crate::fields::traits::{FiniteField, PthRootExtraction};
use crate::polynomials::{DensePolynomial, PolynomialError, traits::UnivariatePolynomial};

impl<F: Field> UnivariatePolynomial<F> for DensePolynomial<F> {
    fn constant(value: F::Elem) -> Self {
        DensePolynomial::constant(value)
    }

    fn degree(&self) -> Option<usize> {
        DensePolynomial::degree(self)
    }

    fn leading_coefficient(&self) -> Option<&F::Elem> {
        DensePolynomial::leading_coefficient(self)
    }

    fn constant_term(&self) -> Option<&F::Elem> {
        DensePolynomial::constant_term(self)
    }

    fn is_zero(&self) -> bool {
        DensePolynomial::is_zero(self)
    }

    fn is_monic(&self) -> bool {
        DensePolynomial::is_monic(self)
    }

    fn add(&self, rhs: &Self) -> Self {
        DensePolynomial::add(self, rhs)
    }

    fn neg(&self) -> Self {
        DensePolynomial::neg(self)
    }

    fn sub(&self, rhs: &Self) -> Self {
        DensePolynomial::sub(self, rhs)
    }

    fn scale(&self, scalar: &F::Elem) -> Self {
        DensePolynomial::scale(self, scalar)
    }

    fn mul(&self, rhs: &Self) -> Self {
        DensePolynomial::mul(self, rhs)
    }

    fn derivative(&self) -> Self {
        DensePolynomial::derivative(self)
    }

    fn gcd(&self, rhs: &Self) -> Self {
        DensePolynomial::gcd(self, rhs)
    }

    fn make_monic(&self) -> Result<Self, PolynomialError> {
        DensePolynomial::make_monic(self)
    }
}

impl<F: FiniteField> PthRootExtraction for DensePolynomial<F>
where
    F::Elem: PthRootExtraction,
{
    fn pth_root(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::new(Vec::new()));
        }

        let characteristic = F::characteristic().to_positive_biguint()?;
        let mut coefficients = Vec::new();

        for (degree, coefficient) in self.coefficients.iter().enumerate() {
            if F::is_zero(coefficient) {
                continue;
            }

            let degree = BigUint::from(degree);
            if &degree % &characteristic != BigUint::ZERO {
                return None;
            }

            let root_degree = (&degree / &characteristic).to_usize()?;
            if coefficients.len() <= root_degree {
                coefficients.resize_with(root_degree + 1, F::zero);
            }

            coefficients[root_degree] = coefficient.pth_root()?;
        }

        Some(Self::new(coefficients))
    }
}

impl<F: FiniteField> DensePolynomial<F> {
    /// Inverts the coordinate substitution `x' ↦ x^p` on one polynomial when
    /// possible.
    ///
    /// This is intentionally different from [`PthRootExtraction`]: the
    /// coefficients are kept fixed and only the exponents are required to lie
    /// in the image of multiplication by `p`.
    ///
    /// Concretely, this returns `Some(P(x'))` exactly when `self` has the form
    /// `P(x^p)`.
    pub(crate) fn inverse_absolute_frobenius_pullback_from_twist(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::new(Vec::new()));
        }

        let characteristic = F::characteristic().to_positive_biguint()?;
        let mut coefficients = Vec::new();

        for (degree, coefficient) in self.coefficients.iter().enumerate() {
            if F::is_zero(coefficient) {
                continue;
            }

            let degree = BigUint::from(degree);
            if &degree % &characteristic != BigUint::ZERO {
                return None;
            }

            let preimage_degree = (&degree / &characteristic).to_usize()?;
            if coefficients.len() <= preimage_degree {
                coefficients.resize_with(preimage_degree + 1, F::zero);
            }

            coefficients[preimage_degree] = coefficient.clone();
        }

        Some(Self::new(coefficients))
    }

    /// Raises one polynomial to a non-negative integer power modulo another.
    ///
    /// This computes `base^exponent mod modulus` by repeated squaring in the
    /// quotient ring `F[x] / (modulus)`.
    ///
    /// The modulus must be non-zero.
    ///
    /// Complexity:
    /// `Θ(log exponent)` polynomial squarings/multiplications, each followed by
    /// one Euclidean remainder reduction modulo `modulus`.
    #[allow(dead_code)]
    pub(crate) fn pow_mod(
        base: &Self,
        exponent: &BigUint,
        modulus: &Self,
    ) -> Result<Self, PolynomialError> {
        if modulus.is_zero() {
            return Err(PolynomialError::DivisionByZeroPolynomial);
        }

        let mut result = Self::constant(F::one());
        let mut factor = base.rem(modulus)?;
        let mut remaining = exponent.clone();

        while !remaining.is_zero() {
            if (&remaining & BigUint::one()) == BigUint::one() {
                result = result.mul(&factor).rem(modulus)?;
            }

            remaining >>= 1usize;
            if !remaining.is_zero() {
                factor = factor.mul(&factor).rem(modulus)?;
            }
        }

        Ok(result)
    }
}
