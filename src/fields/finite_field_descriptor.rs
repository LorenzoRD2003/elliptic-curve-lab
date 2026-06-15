use core::fmt;
use core::num::NonZeroU32;

use crate::fields::error::FieldError;

pub(crate) fn is_valid_field_modulus(modulus: u64) -> bool {
    modulus >= 2
}

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

    /// Returns a compact educational string such as `F_17` or `F_(43^2)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn pretty(&self) -> String {
        let characteristic = self.characteristic;
        let extension_degree = self.extension_degree.get();

        if extension_degree == 1 {
            format!("F_{characteristic}")
        } else {
            format!("F_({characteristic}^{extension_degree})")
        }
    }
}

impl fmt::Display for FiniteFieldDescriptor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.pretty())
    }
}

#[cfg(test)]
mod tests {
    use core::num::NonZeroU32;

    use crate::fields::traits::FiniteField;
    use crate::fields::{FieldError, Fp, finite_field_descriptor::FiniteFieldDescriptor};

    type F17 = Fp<17>;

    #[test]
    fn finite_field_descriptor_pretty_formats_prime_fields() {
        let descriptor = FiniteFieldDescriptor::new(17, NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F17 descriptor should be valid");

        assert_eq!(descriptor.pretty(), "F_17");
        assert_eq!(format!("{descriptor}"), descriptor.pretty());
        assert_eq!(descriptor.cardinality(), Ok(17));
    }

    #[test]
    fn finite_field_descriptor_pretty_formats_extension_fields() {
        let descriptor = FiniteFieldDescriptor::new(43, NonZeroU32::new(2).expect("2 is non-zero"))
            .expect("F43^2 descriptor should be valid");

        assert_eq!(descriptor.pretty(), "F_(43^2)");
        assert_eq!(format!("{descriptor}"), descriptor.pretty());
        assert_eq!(descriptor.cardinality(), Ok(43_u128.pow(2)));
    }

    #[test]
    fn finite_field_descriptor_rejects_invalid_modulus() {
        let error = FiniteFieldDescriptor::new(1, NonZeroU32::new(1).expect("1 is non-zero"))
            .expect_err("characteristic 1 should be rejected");

        assert_eq!(error, FieldError::InvalidModulus { modulus: 1 });
    }

    #[test]
    fn finite_field_trait_descriptor_uses_the_metadata_of_the_type() {
        let descriptor = F17::descriptor().expect("Fp<17> should have a valid descriptor");

        assert_eq!(
            descriptor,
            FiniteFieldDescriptor::new(17, NonZeroU32::new(1).expect("1 is non-zero"))
                .expect("F17 descriptor should be valid")
        );
    }
}
