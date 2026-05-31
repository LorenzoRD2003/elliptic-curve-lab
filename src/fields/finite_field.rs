use core::num::NonZeroU32;

use crate::fields::{errors::FieldError, traits::FiniteField, utils::is_valid_field_modulus};

/// Lightweight metadata describing a finite field family.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteFieldDescriptor {
    /// Prime characteristic of the field.
    pub characteristic: u64,
    /// Extension degree over the prime field.
    ///
    /// The degree is stored as non-zero because a finite field always contains
    /// at least the prime field itself.
    pub extension_degree: NonZeroU32,
}

impl FiniteFieldDescriptor {
    /// Builds a descriptor and performs basic consistency checks.
    pub fn new(characteristic: u64, extension_degree: NonZeroU32) -> Result<Self, FieldError> {
        let descriptor = Self {
            characteristic,
            extension_degree,
        };
        descriptor.check()?;
        Ok(descriptor)
    }

    /// Performs lightweight structural checks only.
    pub fn check(&self) -> Result<(), FieldError> {
        if !is_valid_field_modulus(self.characteristic) {
            return Err(FieldError::InvalidModulus {
                modulus: self.characteristic,
            });
        }

        Ok(())
    }

    /// Computes the field cardinality when it fits into `u128`.
    pub fn cardinality(&self) -> Result<u128, FieldError> {
        u128::from(self.characteristic)
            .checked_pow(self.extension_degree.get())
            .ok_or(FieldError::CardinalityOverflow)
    }
}

/// Creates a descriptor for a concrete finite field type.
pub fn descriptor_for<F: FiniteField>() -> Result<FiniteFieldDescriptor, FieldError> {
    FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
}
