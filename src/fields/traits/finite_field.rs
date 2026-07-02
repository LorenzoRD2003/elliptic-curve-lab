use core::num::NonZeroU32;

use num_bigint::BigUint;

use crate::fields::{
    error::FieldError, finite_field_descriptor::FiniteFieldDescriptor, traits::Field,
};

/// Metadata and validation hooks for finite fields.
pub trait FiniteField: Field {
    /// Returns the degree of the extension over the prime field.
    fn extension_degree() -> NonZeroU32 {
        NonZeroU32::MIN
    }

    /// Returns the field cardinality as an arbitrary-precision integer.
    fn cardinality_biguint() -> BigUint {
        Self::characteristic()
            .to_positive_biguint()
            .expect("finite fields must have positive characteristic")
            .pow(Self::extension_degree().get())
    }

    /// Returns the represented field order `q`.
    ///
    /// This is the checked ergonomic wrapper around
    /// [`Self::cardinality_biguint`]: it first validates the field metadata
    /// through [`Self::check_structure`], then returns the exact field order.
    fn order() -> Result<BigUint, FieldError> {
        Self::check_structure()?;
        Ok(Self::cardinality_biguint())
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

    /// Builds the lightweight finite-field descriptor of this field family.
    fn descriptor() -> Result<FiniteFieldDescriptor, FieldError> {
        let characteristic =
            Self::characteristic()
                .to_positive_biguint()
                .ok_or(FieldError::InvalidModulus {
                    modulus: "0".into(),
                })?;
        FiniteFieldDescriptor::new_biguint(characteristic, Self::extension_degree())
    }
}
