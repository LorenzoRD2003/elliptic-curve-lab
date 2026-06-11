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

/// Runtime-ambient field interface.
///
/// This trait is intentionally parallel to [`Field`], not a replacement for
/// it.
///
/// The existing [`Field`] trait is a static family interface: its identities,
/// arithmetic, characteristic, and small-integer embedding are all determined
/// by the implementing type alone.
///
/// Some mathematically natural fields in this repository depend on additional
/// runtime context. An example is the function field `F(E)` of one concrete
/// short-Weierstrass curve `E`, where arithmetic depends on the chosen cubic
/// relation `y^2 = x^3 + ax + b`.
///
/// For those ambient, runtime-dependent situations, this trait keeps the
/// algebraic operations on the family object itself:
///
/// - identities come from `&self`
/// - arithmetic can validate that two elements are in the same ambient family
/// - division can fail through the family's own error type
pub trait AmbientField {
    /// Stored element type of the ambient field.
    type Elem;

    /// Recoverable failure surface for ambient arithmetic.
    type Error;

    /// Returns the additive identity.
    fn zero(&self) -> Self::Elem;

    /// Returns the multiplicative identity.
    fn one(&self) -> Self::Elem;

    /// Returns whether two elements are equal in the ambient field.
    fn eq(&self, left: &Self::Elem, right: &Self::Elem) -> bool;

    /// Returns whether the given element is zero.
    fn is_zero(&self, value: &Self::Elem) -> bool {
        self.eq(value, &self.zero())
    }

    /// Returns whether the given element is one.
    fn is_one(&self, value: &Self::Elem) -> bool {
        self.eq(value, &self.one())
    }

    /// Adds two elements in the ambient field.
    fn add(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Returns the additive inverse of one element.
    fn neg(&self, value: &Self::Elem) -> Self::Elem;

    /// Subtracts two elements in the ambient field.
    fn sub(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        self.add(left, &self.neg(right))
    }

    /// Multiplies two elements in the ambient field.
    fn mul(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Returns the multiplicative inverse when it exists.
    fn inverse(&self, value: &Self::Elem) -> Result<Self::Elem, Self::Error>;

    /// Divides `left` by `right` when `right` is invertible.
    fn div(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        self.mul(left, &self.inverse(right)?)
    }
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
