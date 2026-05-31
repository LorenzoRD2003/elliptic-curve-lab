use core::num::NonZeroU32;

use crate::fields::errors::FieldError;

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
    /// - [`crate::fields::Q`] returns `0`
    /// - [`crate::fields::ComplexApprox`] returns `0`
    fn characteristic() -> u64;

    fn zero() -> Self::Elem;
    fn one() -> Self::Elem;
    fn from_i64(n: i64) -> Self::Elem;

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
    fn pow(x: &Self::Elem, exponent: u64) -> Self::Elem {
        let mut result = Self::one();
        let mut base = x.clone();
        let mut exp = exponent;

        while exp > 0 {
            if exp & 1 == 1 {
                result = Self::mul(&result, &base);
            }

            exp >>= 1;

            if exp > 0 {
                base = Self::square(&base);
            }
        }

        result
    }

    /// Divides `x` by `y` if `y` is invertible.
    fn div(x: &Self::Elem, y: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Ok(Self::mul(x, &Self::inverse(y)?))
    }

    /// Builds an element from a small integer using the field's canonical map.
    fn elem_from_u64(value: u64) -> Self::Elem;
}

/// Metadata and validation hooks for finite fields.
pub trait FiniteField: Field {
    /// Returns the degree of the extension over the prime field.
    fn extension_degree() -> NonZeroU32 {
        NonZeroU32::MIN
    }

    /// Returns the field cardinality when it fits the chosen representation.
    fn cardinality() -> Option<u128> {
        let characteristic = u128::from(<Self as Field>::characteristic());
        characteristic.checked_pow(Self::extension_degree().get())
    }

    /// Returns whether the field is a prime field.
    fn is_prime_field() -> bool {
        Self::extension_degree().get() == 1
    }

    /// Returns whether the field metadata looks internally consistent.
    fn has_valid_structure() -> bool {
        Self::check_structure().is_ok()
    }

    /// Performs lightweight structural checks for the field family.
    fn check_structure() -> Result<(), FieldError>;

    /// Creates an element from a canonical small integer representation.
    fn try_elem_from_u64(value: u64) -> Result<Self::Elem, FieldError> {
        Self::check_structure()?;
        Ok(Self::elem_from_u64(value))
    }
}

/// Finite fields whose full element set is small enough to enumerate directly.
///
/// This trait is intentionally narrower than [`FiniteField`]. It exists for
/// educational tasks such as exhaustive tables or curve-point enumeration, not
/// as a claim that every finite field backend should always materialize all of
/// its elements.
pub trait EnumerableFiniteField: FiniteField {
    /// Returns every field element in a deterministic order.
    fn elements() -> Vec<Self::Elem>;
}
