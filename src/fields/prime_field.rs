use core::num::NonZeroU32;

use crate::fields::{
    errors::FieldError,
    traits::{Field, FiniteField},
    utils::{extended_gcd, is_prime_u64, is_valid_field_modulus},
};

/// Prime field namespace parameterized by its modulus.
///
/// The associated element type is [`FpElem<P>`], exposed through the
/// [`Field::Elem`] associated type.
#[derive(Clone, Copy, Debug)]
pub struct Fp<const P: u64>;

/// Canonical residue class modulo `P`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct FpElem<const P: u64> {
    value: u64,
}

impl<const P: u64> Fp<P> {
    /// Returns the modulus of the field.
    pub const fn modulus() -> u64 {
        P
    }

    /// Performs structural checks for the prime field.
    ///
    /// This is the central validation hook for the field family. Arithmetic
    /// operations intentionally assume the modulus was already validated through
    /// construction or explicit structure checks.
    pub fn validate_modulus() -> Result<(), FieldError> {
        if !is_valid_field_modulus(P) {
            return Err(FieldError::InvalidModulus { modulus: P });
        }

        if !is_prime_u64(P) {
            return Err(FieldError::InvalidModulus { modulus: P });
        }

        Ok(())
    }

    /// Builds an element by reducing the input modulo `P`.
    ///
    /// This is a small convenience wrapper around [`FpElem::new`].
    pub fn new_elem(value: u64) -> Result<FpElem<P>, FieldError> {
        FpElem::new(value)
    }

    /// Reduces an unsigned integer modulo `P`.
    ///
    /// This helper assumes the modulus was already validated by construction or
    /// explicit structure checks.
    fn reduce_u128(value: u128) -> u64 {
        (value % u128::from(P)) as u64
    }

    /// Reduces a signed integer modulo `P` using Euclidean remainder semantics.
    ///
    /// This helper is primarily used by [`Field::from_i64`].
    fn reduce_i64(value: i64) -> u64 {
        let modulus = i128::from(P.max(1));
        let reduced = i128::from(value).rem_euclid(modulus);
        reduced as u64
    }

    /// Computes the multiplicative inverse of a non-zero canonical residue.
    ///
    /// The returned value is also canonical modulo `P`.
    fn inverse_value(value: u64) -> Option<u64> {
        let (gcd, coefficient, _) = extended_gcd(i128::from(value), i128::from(P));
        if gcd != 1 {
            return None;
        }

        let inverse = coefficient.rem_euclid(i128::from(P)) as u64;
        Some(inverse)
    }
}

impl<const P: u64> FpElem<P> {
    /// Builds an element and reduces the input into canonical form.
    pub fn new(value: u64) -> Result<Self, FieldError> {
        Fp::<P>::validate_modulus()?;
        Ok(Self {
            value: Fp::<P>::reduce_u128(u128::from(value)),
        })
    }

    /// Builds an element only if the value is already canonical.
    ///
    /// This constructor validates the modulus and rejects representatives that
    /// are not already in the interval `[0, P)`.
    pub fn try_from_canonical(value: u64) -> Result<Self, FieldError> {
        Fp::<P>::validate_modulus()?;
        if value >= P {
            return Err(FieldError::ElementOutOfRange {
                value: value.to_string(),
            });
        }

        Ok(Self { value })
    }

    /// Returns the modulus associated with this element type.
    pub const fn modulus() -> u64 {
        P
    }

    /// Returns the canonical representative.
    pub const fn value(&self) -> u64 {
        self.value
    }

    /// Returns whether the stored representative is canonical.
    pub fn is_canonical(&self) -> bool {
        self.value < P
    }
}

impl<const P: u64> Field for Fp<P> {
    /// Finite prime fields are not algebraically closed.
    ///
    /// They admit non-trivial algebraic closures, but `Fp<P>` models only the
    /// prime field itself, not that larger closure.
    const IS_ALGEBRAICALLY_CLOSED: bool = false;

    type Elem = FpElem<P>;

    /// Returns the additive identity.
    fn zero() -> Self::Elem {
        Self::Elem { value: 0 }
    }

    /// Returns the multiplicative identity.
    fn one() -> Self::Elem {
        if P >= 2 {
            Self::Elem { value: 1 % P }
        } else {
            Self::Elem { value: 0 }
        }
    }

    /// Embeds a signed integer into the field using the canonical residue map.
    fn from_i64(n: i64) -> Self::Elem {
        Self::Elem {
            value: Self::reduce_i64(n),
        }
    }

    /// Adds two canonical residues modulo `P`.
    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::Elem {
            value: Self::reduce_u128(u128::from(x.value) + u128::from(y.value)),
        }
    }

    /// Subtracts `y` from `x` modulo `P`.
    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::Elem {
            value: Self::reduce_u128(u128::from(x.value) + u128::from(P) - u128::from(y.value)),
        }
    }

    /// Multiplies two canonical residues modulo `P`.
    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        FpElem {
            value: Fp::<P>::reduce_u128(u128::from(x.value) * u128::from(y.value)),
        }
    }

    /// Returns the additive inverse of `x`.
    fn neg(x: &Self::Elem) -> Self::Elem {
        if x.value == 0 {
            return *x;
        }

        Self::Elem { value: P - x.value }
    }

    /// Returns the multiplicative inverse of `x` when `x` is non-zero.
    ///
    /// This method assumes the field modulus is already valid.
    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        if x.value == 0 {
            return None;
        }

        Self::inverse_value(x.value).map(|value| FpElem { value })
    }

    /// Returns whether the canonical representatives are equal.
    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        x.value == y.value
    }

    /// Returns the multiplicative inverse of `x` or a structured error.
    ///
    /// This method does not revalidate the modulus on every call; it relies on
    /// the field type having been validated when values were constructed.
    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        if x.value == 0 {
            return Err(FieldError::DivisionByZero);
        }

        Self::inv(x).ok_or(FieldError::Unsupported(
            "non-zero prime-field element should always be invertible",
        ))
    }

    /// Embeds an unsigned integer into the field using modular reduction.
    ///
    /// This method assumes the field modulus was validated elsewhere.
    fn elem_from_u64(value: u64) -> Self::Elem {
        Self::Elem {
            value: Self::reduce_u128(u128::from(value)),
        }
    }
}

impl<const P: u64> FiniteField for Fp<P> {
    /// Returns the characteristic of the field.
    fn characteristic() -> u64 {
        P
    }

    /// Returns the extension degree over the prime field.
    fn extension_degree() -> NonZeroU32 {
        NonZeroU32::MIN
    }

    /// Verifies that `P` can define a prime field.
    fn check_structure() -> Result<(), FieldError> {
        Self::validate_modulus()
    }
}

#[cfg(test)]
mod tests {
    use std::hint::black_box;

    use super::{Fp, FpElem};
    use crate::fields::{Field, FieldError, FiniteField};

    type F17 = Fp<17>;
    type E17 = FpElem<17>;

    fn e(value: u64) -> E17 {
        F17::elem_from_u64(value)
    }

    fn all_elements() -> Vec<E17> {
        (0..17).map(e).collect()
    }

    #[test]
    fn constructor_reduces_values() {
        assert_eq!(E17::new(0).expect("valid element").value(), 0);
        assert_eq!(E17::new(17).expect("valid element").value(), 0);
        assert_eq!(E17::new(18).expect("valid element").value(), 1);
        assert_eq!(E17::new(52).expect("valid element").value(), 1);
    }

    #[test]
    fn canonical_constructor_accepts_and_rejects_as_expected() {
        assert_eq!(
            E17::try_from_canonical(16)
                .expect("canonical value should be accepted")
                .value(),
            16
        );

        let error = E17::try_from_canonical(17).expect_err("non-canonical value should fail");
        assert!(matches!(error, FieldError::ElementOutOfRange { .. }));
    }

    #[test]
    fn zero_one_and_modulus_are_exposed() {
        assert_eq!(F17::modulus(), 17);
        assert_eq!(E17::modulus(), 17);
        assert_eq!(F17::zero().value(), 0);
        assert_eq!(F17::one().value(), 1);
        assert!(F17::is_zero(&F17::zero()));
        assert!(!F17::is_zero(&F17::one()));
    }

    #[test]
    fn from_i64_handles_negative_inputs() {
        assert_eq!(F17::from_i64(0).value(), 0);
        assert_eq!(F17::from_i64(19).value(), 2);
        assert_eq!(F17::from_i64(-1).value(), 16);
        assert_eq!(F17::from_i64(-18).value(), 16);
        assert_eq!(F17::from_i64(-34).value(), 0);
    }

    #[test]
    fn addition_subtraction_and_negation_wrap_correctly() {
        assert_eq!(F17::add(&e(9), &e(11)).value(), 3);
        assert_eq!(F17::sub(&e(3), &e(5)).value(), 15);
        assert_eq!(F17::sub(&e(5), &e(3)).value(), 2);
        assert_eq!(F17::neg(&e(0)).value(), 0);
        assert_eq!(F17::neg(&e(1)).value(), 16);
        assert_eq!(F17::neg(&e(9)).value(), 8);
    }

    #[test]
    fn multiplication_square_and_cube_work() {
        assert_eq!(F17::mul(&e(5), &e(7)).value(), 1);
        assert_eq!(F17::square(&e(6)).value(), 2);
        assert_eq!(F17::cube(&e(3)).value(), 10);
    }

    #[test]
    fn pow_handles_edge_cases_and_fermat_relations() {
        assert_eq!(F17::pow(&e(5), 0).value(), 1);
        assert_eq!(F17::pow(&e(5), 1).value(), 5);
        assert_eq!(F17::pow(&e(5), 2).value(), 8);
        assert_eq!(F17::pow(&e(5), 16).value(), 1);
        assert_eq!(F17::pow(&e(5), 17).value(), 5);
        assert_eq!(F17::pow(&e(0), 7).value(), 0);
    }

    #[test]
    fn option_inverse_and_result_inverse_agree() {
        assert_eq!(F17::inv(&e(0)), None);
        assert!(matches!(
            F17::inverse(&e(0)),
            Err(FieldError::DivisionByZero)
        ));

        for value in 1..17 {
            let x = e(value);
            let inv_option = F17::inv(&x).expect("non-zero element should be invertible");
            let inv_result = F17::inverse(&x).expect("non-zero element should be invertible");
            assert_eq!(inv_option, inv_result);
            assert_eq!(F17::mul(&x, &inv_result), F17::one());
        }
    }

    #[test]
    fn division_matches_multiplication_by_inverse() {
        let quotient = F17::div(&e(8), &e(3)).expect("division should succeed");
        let expected = F17::mul(&e(8), &F17::inverse(&e(3)).expect("inverse should exist"));

        assert_eq!(quotient, expected);
        assert_eq!(F17::mul(&quotient, &e(3)), e(8));
        assert!(matches!(
            F17::div(&e(8), &e(0)),
            Err(FieldError::DivisionByZero)
        ));
    }

    #[test]
    fn equality_and_zero_checks_match_canonical_values() {
        assert!(<F17 as Field>::eq(&e(2), &e(19)));
        assert!(!<F17 as Field>::eq(&e(2), &e(3)));
        assert!(F17::is_zero(&e(17)));
        assert!(!F17::is_zero(&e(16)));
    }

    #[test]
    fn finite_field_metadata_is_correct() {
        assert_eq!(F17::characteristic(), 17);
        assert_eq!(F17::extension_degree().get(), 1);
        assert_eq!(F17::cardinality(), Some(17));
        assert!(F17::is_prime_field());
        assert!(F17::has_valid_structure());
        assert!(!black_box(F17::IS_ALGEBRAICALLY_CLOSED));
        assert_eq!(F17::try_elem_from_u64(20).expect("field is valid"), e(20));
    }

    #[test]
    fn invalid_moduli_are_rejected() {
        assert!(matches!(
            Fp::<0>::check_structure(),
            Err(FieldError::InvalidModulus { modulus: 0 })
        ));
        assert!(matches!(
            Fp::<1>::check_structure(),
            Err(FieldError::InvalidModulus { modulus: 1 })
        ));
        assert!(matches!(
            Fp::<15>::check_structure(),
            Err(FieldError::InvalidModulus { modulus: 15 })
        ));
    }

    #[test]
    fn exhaustive_additive_group_properties_hold_for_f17() {
        let elements = all_elements();

        for x in &elements {
            assert_eq!(F17::add(x, &F17::zero()), *x);
            assert_eq!(F17::add(&F17::zero(), x), *x);
            assert_eq!(F17::add(x, &F17::neg(x)), F17::zero());
        }

        for x in &elements {
            for y in &elements {
                assert_eq!(F17::add(x, y), F17::add(y, x));
                assert_eq!(F17::sub(x, y), F17::add(x, &F17::neg(y)));
            }
        }

        for x in &elements {
            for y in &elements {
                for z in &elements {
                    let left = F17::add(&F17::add(x, y), z);
                    let right = F17::add(x, &F17::add(y, z));
                    assert_eq!(left, right);
                }
            }
        }
    }

    #[test]
    fn exhaustive_multiplicative_properties_hold_for_f17() {
        let elements = all_elements();

        for x in &elements {
            assert_eq!(F17::mul(x, &F17::one()), *x);
            assert_eq!(F17::mul(&F17::one(), x), *x);
        }

        for x in &elements {
            for y in &elements {
                assert_eq!(F17::mul(x, y), F17::mul(y, x));
            }
        }

        for x in &elements {
            for y in &elements {
                for z in &elements {
                    let left = F17::mul(&F17::mul(x, y), z);
                    let right = F17::mul(x, &F17::mul(y, z));
                    assert_eq!(left, right);
                }
            }
        }

        for value in 1..17 {
            let x = e(value);
            let inverse = F17::inverse(&x).expect("non-zero element should be invertible");
            assert_eq!(F17::mul(&x, &inverse), F17::one());
            assert_eq!(F17::mul(&inverse, &x), F17::one());
        }
    }

    #[test]
    fn exhaustive_distributivity_holds_for_f17() {
        let elements = all_elements();

        for x in &elements {
            for y in &elements {
                for z in &elements {
                    let left = F17::mul(x, &F17::add(y, z));
                    let right = F17::add(&F17::mul(x, y), &F17::mul(x, z));
                    assert_eq!(left, right);
                }
            }
        }
    }
}
