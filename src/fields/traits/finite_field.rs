use core::num::NonZeroU32;

use crate::fields::{
    error::FieldError, finite_field_descriptor::FiniteFieldDescriptor, traits::Field,
};

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

    /// Builds the lightweight finite-field descriptor of this field family.
    fn descriptor() -> Result<FiniteFieldDescriptor, FieldError> {
        FiniteFieldDescriptor::new(Self::characteristic(), Self::extension_degree())
    }
}
