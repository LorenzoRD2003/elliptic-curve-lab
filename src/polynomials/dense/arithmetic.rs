use crate::fields::traits::*;
use num_bigint::BigInt;

use super::DensePolynomial;

impl<F: Field> DensePolynomial<F> {
    /// Adds two dense polynomials coefficient-wise.
    ///
    /// The result length is the maximum of the two input lengths. Coefficients
    /// beyond the shorter input are treated as zero. The result is normalized
    /// through the constructor, so trailing zeros are trimmed away.
    pub fn add(&self, rhs: &Self) -> Self {
        let target_len = self.len().max(rhs.len());
        let mut coefficients = Vec::with_capacity(target_len);

        for index in 0..target_len {
            let lhs_coeff = self
                .coefficients
                .get(index)
                .cloned()
                .unwrap_or_else(F::zero);
            let rhs_coeff = rhs.coefficients.get(index).cloned().unwrap_or_else(F::zero);
            coefficients.push(F::add(&lhs_coeff, &rhs_coeff));
        }

        Self::new(coefficients)
    }

    /// Negates every coefficient of the polynomial.
    pub fn neg(&self) -> Self {
        Self::new(self.coefficients.iter().map(F::neg).collect())
    }

    /// Subtracts two dense polynomials coefficient-wise.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// Multiplies every coefficient of the polynomial by the same field
    /// element.
    pub fn scale(&self, scalar: &F::Elem) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .map(|coefficient| F::mul(coefficient, scalar))
                .collect(),
        )
    }

    /// Returns the formal derivative of the polynomial. If
    /// `f(x) = a0 + a1*x + a2*x^2 + ... + an*x^n`, then this
    ///  method returns `f'(x) = a1 + 2*a2*x + ... + n*an*x^(n-1)`.
    pub fn derivative(&self) -> Self {
        if self.len() <= 1 {
            return Self::new(Vec::new());
        }

        let coefficients = self
            .coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(degree, coefficient)| {
                let scalar = F::from_bigint(&BigInt::from(degree));
                F::mul(coefficient, &scalar)
            })
            .collect();

        Self::new(coefficients)
    }

    /// Multiplies two dense polynomials using the naive quadratic algorithm.
    ///
    /// This implementation is intentionally straightforward and educational,
    /// and is not optimized.
    ///
    /// If either side is represented by an empty coefficient vector, the
    /// result keeps that empty zero representation.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::new(Vec::new());
        }

        let mut coefficients = vec![F::zero(); self.len() + rhs.len() - 1];

        for (lhs_degree, lhs_coeff) in self.coefficients.iter().enumerate() {
            for (rhs_degree, rhs_coeff) in rhs.coefficients.iter().enumerate() {
                let index = lhs_degree + rhs_degree;
                let term = F::mul(lhs_coeff, rhs_coeff);
                let updated = F::add(&coefficients[index], &term);
                coefficients[index] = updated;
            }
        }

        Self::new(coefficients)
    }

    pub fn square(&self) -> Self {
        self.mul(self)
    }

    pub fn cube(&self) -> Self {
        self.square().mul(self)
    }
}
