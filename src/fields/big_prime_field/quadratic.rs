use num_bigint::BigUint;
use num_traits::One;

use super::BigPrimeField;
use crate::fields::{
    FieldError,
    big_prime_field::BigPrimeFieldElem,
    traits::{AmbientField, QuadraticCharacterValue},
};

impl BigPrimeField {
    /// Evaluates the quadratic character `χ : F_p → {0, ±1}` for odd prime `p`.
    ///
    /// The current implementation follows Euler's criterion:
    ///
    /// `χ(x) = x^((p-1)/2)` for `x ≠ 0`.
    ///
    /// Characteristic `2` is reported as unsupported, matching the existing
    /// trait-side finite-field story elsewhere in the repo.
    pub fn quadratic_character_of(
        &self,
        value: &BigPrimeFieldElem,
    ) -> Result<QuadraticCharacterValue, FieldError> {
        if self.is_zero(value) {
            return Ok(QuadraticCharacterValue::Zero);
        }

        if self.modulus() == &BigUint::from(2u8) {
            return Err(FieldError::Unsupported(
                "quadratic character is only implemented for finite fields of odd characteristic",
            ));
        }

        let exponent = (self.modulus() - BigUint::one()) >> 1usize;
        let character_power = self.pow_biguint(value, &exponent);

        if AmbientField::eq(self, &character_power, &self.one()) {
            return Ok(QuadraticCharacterValue::Residue);
        }

        let minus_one = self.neg(&self.one());
        if AmbientField::eq(self, &character_power, &minus_one) {
            return Ok(QuadraticCharacterValue::NonResidue);
        }

        Err(FieldError::Unsupported(
            "quadratic-character power did not land in {0, ±1} as expected over the current big prime-field backend",
        ))
    }

    /// Returns whether `value` is a quadratic residue in `F_p`.
    pub fn has_square_root(&self, value: &BigPrimeFieldElem) -> Result<bool, FieldError> {
        Ok(matches!(
            self.quadratic_character_of(value)?,
            QuadraticCharacterValue::Zero | QuadraticCharacterValue::Residue
        ))
    }
}
