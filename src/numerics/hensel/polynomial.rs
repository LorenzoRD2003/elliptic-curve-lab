use num_bigint::BigInt;
use num_traits::Zero;

use crate::polynomials::IntegerPolynomial;

pub(in crate::numerics::hensel) trait HenselPolynomialEvaluator {
    fn dense_coefficients(&self) -> Vec<BigInt>;

    fn evaluate(&self, x: &BigInt) -> BigInt;

    fn evaluate_derivative(&self, x: &BigInt) -> BigInt;

    fn is_zero_polynomial(&self) -> bool;

    fn is_constant_polynomial(&self) -> bool;
}

impl HenselPolynomialEvaluator for [BigInt] {
    fn dense_coefficients(&self) -> Vec<BigInt> {
        self.to_vec()
    }

    fn evaluate(&self, x: &BigInt) -> BigInt {
        evaluate_polynomial(self, x)
    }

    fn evaluate_derivative(&self, x: &BigInt) -> BigInt {
        evaluate_derivative(self, x)
    }

    fn is_zero_polynomial(&self) -> bool {
        self.is_empty()
    }

    fn is_constant_polynomial(&self) -> bool {
        self.len() == 1
    }
}

impl<const N: usize> HenselPolynomialEvaluator for [BigInt; N] {
    fn dense_coefficients(&self) -> Vec<BigInt> {
        self.as_slice().dense_coefficients()
    }

    fn evaluate(&self, x: &BigInt) -> BigInt {
        HenselPolynomialEvaluator::evaluate(self.as_slice(), x)
    }

    fn evaluate_derivative(&self, x: &BigInt) -> BigInt {
        HenselPolynomialEvaluator::evaluate_derivative(self.as_slice(), x)
    }

    fn is_zero_polynomial(&self) -> bool {
        self.as_slice().is_zero_polynomial()
    }

    fn is_constant_polynomial(&self) -> bool {
        self.as_slice().is_constant_polynomial()
    }
}

impl HenselPolynomialEvaluator for IntegerPolynomial {
    fn dense_coefficients(&self) -> Vec<BigInt> {
        self.to_dense_coefficients()
    }

    fn evaluate(&self, x: &BigInt) -> BigInt {
        self.evaluate(x)
    }

    fn evaluate_derivative(&self, x: &BigInt) -> BigInt {
        self.evaluate_derivative(x)
    }

    fn is_zero_polynomial(&self) -> bool {
        self.is_zero()
    }

    fn is_constant_polynomial(&self) -> bool {
        self.is_constant()
    }
}

pub(super) fn evaluate_polynomial(coefs: &[BigInt], x: &BigInt) -> BigInt {
    coefs
        .iter()
        .rev()
        .fold(BigInt::zero(), |acc, c| acc * x + c)
}

pub(super) fn evaluate_derivative(coefs: &[BigInt], x: &BigInt) -> BigInt {
    coefs
        .iter()
        .enumerate()
        .skip(1)
        .rev()
        .fold(BigInt::zero(), |acc, (deg, c)| {
            acc * x + c * BigInt::from(deg)
        })
}
