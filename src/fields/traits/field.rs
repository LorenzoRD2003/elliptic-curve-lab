use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

use crate::fields::{FieldCharacteristic, FieldError};

/// The trait stays explicit instead of overloading the standard arithmetic
/// operators so implementors can keep error handling and canonicalization
/// visible while the project is still stabilizing its algebraic APIs.
pub trait Field {
    /// Whether the modeled field is algebraically closed. This is a mathematical
    /// property of the field family represented by the backend.
    ///
    /// Intuition:
    ///
    /// - an algebraically closed field contains a root for every non-constant
    ///   polynomial with coefficients in that field
    /// - a field that is not algebraically closed can still admit useful
    ///   algebraic extensions, and that distinction will matter for later
    ///   APIs around irreducibility and extension construction
    ///
    /// Examples in this crate:
    ///
    /// - [`crate::fields::ComplexApprox`] models `C` approximately, so it sets
    ///   this to `true`
    /// - [`crate::fields::Q`] models the rationals, which are not algebraically
    ///   closed, so it sets this to `false`
    /// - [`crate::fields::Fp`] models finite prime fields, which are not
    ///   algebraically closed, so it sets this to `false`
    const IS_ALGEBRAICALLY_CLOSED: bool;

    type Elem: Clone + std::fmt::Debug;

    /// Returns the characteristic of the field family.
    ///
    /// Examples in this crate:
    ///
    /// - [`crate::fields::Fp`] returns its prime modulus
    /// - [`crate::fields::Q`] returns characteristic zero
    /// - [`crate::fields::ComplexApprox`] returns characteristic zero
    fn characteristic() -> FieldCharacteristic;

    /// Returns whether the field has the given characteristic.
    fn has_characteristic(value: u8) -> bool {
        let value = BigUint::from(value);
        match Self::characteristic() {
            FieldCharacteristic::Zero => value.is_zero(),
            FieldCharacteristic::Positive(characteristic) => characteristic == value,
        }
    }

    fn zero() -> Self::Elem;
    fn one() -> Self::Elem;
    /// Embeds an integer into the field through the canonical map `ℤ → F`.
    fn from_bigint(n: &BigInt) -> Self::Elem;

    /// Embeds a small signed integer into the field.
    ///
    /// This is a literal-convenience wrapper; generic arithmetic should prefer
    /// [`Self::from_bigint`] so it does not bake machine integer widths into
    /// field APIs.
    fn from_i64(n: i64) -> Self::Elem {
        Self::from_bigint(&BigInt::from(n))
    }

    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn neg(x: &Self::Elem) -> Self::Elem;
    fn inv(x: &Self::Elem) -> Option<Self::Elem>;
    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool;

    fn is_zero(x: &Self::Elem) -> bool {
        Self::eq(x, &Self::zero())
    }

    fn square(x: &Self::Elem) -> Self::Elem {
        Self::mul(x, x)
    }

    fn cube(x: &Self::Elem) -> Self::Elem {
        let x2 = Self::square(x);
        Self::mul(&x2, x)
    }

    /// Returns the multiplicative inverse when it exists.
    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError>;

    /// Raises the element to an unsigned power.
    fn pow(x: &Self::Elem, exponent: &BigUint) -> Self::Elem {
        let mut result = Self::one();
        let mut base = x.clone();
        let mut exp = exponent.clone();

        while !exp.is_zero() {
            if (&exp & BigUint::one()) == BigUint::one() {
                result = Self::mul(&result, &base);
            }
            exp >>= 1usize;

            if !exp.is_zero() {
                base = Self::square(&base);
            }
        }
        result
    }

    /// Divides `x` by `y` if `y` is invertible.
    fn div(x: &Self::Elem, y: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Ok(Self::mul(x, &Self::inverse(y)?))
    }
}
