use num_bigint::BigUint;

use crate::fields::{
    ExtensionField, ExtensionFieldSpec, FieldError, FiniteField, Fp,
    pth_root_extraction::finite_field_pow_biguint,
};

/// For a finite field `F_q` of odd characteristic, the quadratic character
/// `χ : F_q -> {0, ±1}` is defined by
///
/// - `χ(0) = 0`
/// - `χ(x) = 1` when `x` is a non-zero square
/// - `χ(x) = -1` when `x` is a non-square
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QuadraticCharacterValue {
    Zero,
    Residue,
    NonResidue,
}

impl QuadraticCharacterValue {
    /// Returns the signed value in `{0, 1, -1}`.
    pub fn as_i8(self) -> i8 {
        match self {
            Self::Zero => 0,
            Self::Residue => 1,
            Self::NonResidue => -1,
        }
    }

    /// Returns the signed value in `{0, 1, -1}` as `i128`.
    pub fn as_i128(self) -> i128 {
        i128::from(self.as_i8())
    }
}

/// Finite fields that can evaluate the quadratic character `χ`.
///
/// The current implementation is intentionally honest about scope:
///
/// - it supports finite fields of odd characteristic
/// - it uses the finite-field identity
///   `χ(x) = x^((q-1)/2)` for `x != 0`
/// - it reports characteristic-`2` backends as unsupported, since the usual
///   three-valued quadratic-character story collapses there
///
/// This trait is field-family-oriented rather than value-oriented because the
/// meaning of `χ` depends on the ambient finite field `F_q`.
pub trait QuadraticCharacterFiniteField: FiniteField + Sized {
    /// Evaluates the quadratic character in a finite field of odd characteristic.
    ///
    /// The current implementation uses the Euler-style finite-field identity
    /// `χ(x) = x^((q-1)/2)` for `x != 0`.
    ///
    /// For even characteristic, it returns a [FieldError].
    ///
    /// Complexity: `Θ(log q)` field multiplications and squarings under the
    /// current repeated-squaring backend.
    fn quadratic_character_of(x: &Self::Elem) -> Result<QuadraticCharacterValue, FieldError> {
        if Self::is_zero(x) {
            return Ok(QuadraticCharacterValue::Zero);
        }

        if Self::characteristic() == 2 {
            return Err(FieldError::Unsupported(
                "quadratic character is only implemented for finite fields of odd characteristic",
            ));
        }

        let field_order = Self::cardinality().ok_or(FieldError::CardinalityOverflow)?;
        let exponent = BigUint::from((field_order - 1) / 2);
        let value = finite_field_pow_biguint::<Self>(x, &exponent);

        if Self::eq(&value, &Self::one()) {
            return Ok(QuadraticCharacterValue::Residue);
        }

        let minus_one = Self::neg(&Self::one());
        if Self::eq(&value, &minus_one) {
            return Ok(QuadraticCharacterValue::NonResidue);
        }

        Err(FieldError::Unsupported(
            "quadratic-character power did not land in {0, ±1} as expected over the current finite-field backend",
        ))
    }
}

impl<const P: u64> QuadraticCharacterFiniteField for Fp<P> {}

impl<S: ExtensionFieldSpec> QuadraticCharacterFiniteField for ExtensionField<S> where
    ExtensionField<S>: FiniteField
{
}

#[cfg(test)]
mod tests {
    use crate::fields::{
        EnumerableFiniteField, Field, Fp, QuadraticCharacterFiniteField, QuadraticCharacterValue,
        SqrtField,
    };

    type F17 = Fp<17>;
    type F2 = Fp<2>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F17Sqrt3QuadraticCharacterSpec,
        field: F17Sqrt3QuadraticCharacter,
        base: F17,
        non_residue: 3,
        name: "F17(sqrt(3)) for quadratic-character tests",
    );

    #[test]
    fn prime_field_quadratic_character_distinguishes_zero_residue_and_non_residue() {
        assert_eq!(
            F17::quadratic_character_of(&F17::zero()),
            Ok(QuadraticCharacterValue::Zero)
        );
        assert_eq!(
            F17::quadratic_character_of(&F17::from_i64(4)),
            Ok(QuadraticCharacterValue::Residue)
        );
        assert_eq!(
            F17::quadratic_character_of(&F17::from_i64(3)),
            Ok(QuadraticCharacterValue::NonResidue)
        );
    }

    #[test]
    fn prime_field_quadratic_character_matches_square_root_existence_away_from_zero() {
        for element in F17::elements() {
            let character = F17::quadratic_character_of(&element)
                .expect("odd prime fields should support the quadratic character");
            if F17::is_zero(&element) {
                assert_eq!(character, QuadraticCharacterValue::Zero);
            } else if F17::has_square_root(&element) {
                assert_eq!(character, QuadraticCharacterValue::Residue);
            } else {
                assert_eq!(character, QuadraticCharacterValue::NonResidue);
            }
        }
    }

    #[test]
    fn extension_field_quadratic_character_matches_square_root_existence_away_from_zero() {
        for element in F17Sqrt3QuadraticCharacter::elements() {
            let character = F17Sqrt3QuadraticCharacter::quadratic_character_of(&element)
                .expect("odd finite extensions should support the quadratic character");
            if F17Sqrt3QuadraticCharacter::is_zero(&element) {
                assert_eq!(character, QuadraticCharacterValue::Zero);
            } else if F17Sqrt3QuadraticCharacter::has_square_root(&element) {
                assert_eq!(character, QuadraticCharacterValue::Residue);
            } else {
                assert_eq!(character, QuadraticCharacterValue::NonResidue);
            }
        }
    }

    #[test]
    fn characteristic_two_prime_field_is_reported_as_unsupported() {
        assert_eq!(
            F2::quadratic_character_of(&F2::one()),
            Err(crate::fields::FieldError::Unsupported(
                "quadratic character is only implemented for finite fields of odd characteristic"
            ))
        );
    }

    #[test]
    fn quadratic_character_value_converts_to_signed_integers() {
        assert_eq!(QuadraticCharacterValue::Zero.as_i8(), 0);
        assert_eq!(QuadraticCharacterValue::Residue.as_i8(), 1);
        assert_eq!(QuadraticCharacterValue::NonResidue.as_i8(), -1);
        assert_eq!(QuadraticCharacterValue::NonResidue.as_i128(), -1);
    }
}
