use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::BigPrimeField;
use crate::fields::{big_prime_field::BigPrimeFieldElem, traits::AmbientField};

impl BigPrimeField {
    /// Raises one residue to a non-negative `BigUint` exponent in `F_p`.
    ///
    /// This follows the usual repeated-squaring story in the ambient field:
    ///
    /// - start from `1`
    /// - square the running base each step
    /// - multiply into the accumulator when the current exponent bit is `1`
    ///
    /// Complexity: `Θ(log n)` field operations for exponent `n`.
    pub fn pow_biguint(&self, value: &BigPrimeFieldElem, exponent: &BigUint) -> BigPrimeFieldElem {
        let mut result = self.one();
        let mut base = value.clone();
        let mut remaining_exponent = exponent.clone();

        while !remaining_exponent.is_zero() {
            if (&remaining_exponent & BigUint::one()) == BigUint::one() {
                result = self
                    .mul(&result, &base)
                    .expect("prime-field multiplication should stay total");
            }

            remaining_exponent >>= 1usize;

            if !remaining_exponent.is_zero() {
                base = self
                    .mul(&base, &base)
                    .expect("prime-field multiplication should stay total");
            }
        }
        result
    }
}
