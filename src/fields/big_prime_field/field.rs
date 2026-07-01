use core::fmt;

use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;
use num_traits::{One, Zero};

use crate::fields::{FieldError, big_prime_field::BigPrimeFieldElem, traits::AmbientField};
use crate::numerics::inverse_mod_biguint;

/// Runtime prime field `F_p` for primes that do not naturally fit the current
/// type-level `Fp<P>` surface.
///
/// This runtime-owned family is the large-prime analogue of the existing
/// static `Fp<P>` backend:
///
/// - the prime modulus `p` is stored as a `BigUint`
/// - elements are canonical residues in `[0, p)`
/// - arithmetic is performed relative to `&self`
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BigPrimeField {
    modulus: BigUint,
}

impl fmt::Display for BigPrimeField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "F_{}", self.modulus)
    }
}

impl BigPrimeField {
    /// Builds one runtime prime field `F_p`.
    ///
    /// The current validation is intentionally honest rather than clever:
    /// `p` must be at least `2` and must pass `num-prime`'s primality check.
    pub fn new(modulus: BigUint) -> Result<Self, FieldError> {
        validate_prime_modulus(&modulus)?;
        Ok(Self { modulus })
    }

    /// Returns the prime modulus `p`.
    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }

    /// Returns the bit-length of the prime modulus `p`.
    pub fn modulus_bits(&self) -> u64 {
        self.modulus.bits()
    }

    /// Reduces one non-negative integer into its canonical residue class.
    pub fn elem(&self, value: BigUint) -> BigPrimeFieldElem {
        BigPrimeFieldElem::new_canonical(value % &self.modulus)
    }

    /// Builds a canonical residue class from a small unsigned integer.
    pub fn elem_from_u64(&self, value: u64) -> BigPrimeFieldElem {
        self.elem(BigUint::from(value))
    }

    /// Returns whether `value` is already canonical for this field.
    pub fn is_canonical_value(&self, value: &BigUint) -> bool {
        value < &self.modulus
    }
}

impl AmbientField for BigPrimeField {
    type Elem = BigPrimeFieldElem;
    type Error = FieldError;

    fn zero(&self) -> Self::Elem {
        BigPrimeFieldElem::new_canonical(BigUint::zero())
    }

    fn one(&self) -> Self::Elem {
        BigPrimeFieldElem::new_canonical(BigUint::one())
    }

    fn eq(&self, left: &Self::Elem, right: &Self::Elem) -> bool {
        left == right
    }

    fn add(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        let sum_in_z = left.value() + right.value();
        Ok(self.elem(sum_in_z))
    }

    fn neg(&self, value: &Self::Elem) -> Self::Elem {
        if self.is_zero(value) {
            return self.zero();
        }
        let additive_inverse_in_z = self.modulus() - value.value();
        BigPrimeFieldElem::new_canonical(additive_inverse_in_z)
    }

    fn mul(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        let product_in_z = left.value() * right.value();
        Ok(self.elem(product_in_z))
    }

    fn inverse(&self, value: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        if self.is_zero(value) {
            return Err(FieldError::DivisionByZero);
        }
        let inverse_mod_p = inverse_mod_biguint(value.value(), &self.modulus)
            .expect("non-zero residues in a prime field should be invertible");
        Ok(BigPrimeFieldElem::new_canonical(inverse_mod_p))
    }
}

fn validate_prime_modulus(modulus: &BigUint) -> Result<(), FieldError> {
    if modulus < &BigUint::from(2u8) || !is_prime(modulus, None).probably() {
        return Err(FieldError::InvalidBigModulus {
            modulus: modulus.to_string(),
        });
    }
    Ok(())
}
