use core::{fmt, num::NonZeroU32};

use num_bigint::{BigInt, BigUint};

use crate::fields::error::FieldError;

fn is_valid_big_field_modulus(modulus: &BigUint) -> bool {
    modulus >= &BigUint::from(2u8)
}

/// Lightweight metadata describing a finite field family.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteFieldDescriptor {
    /// Prime characteristic of the field.
    pub characteristic: BigUint,
    /// Extension degree over the prime field.
    ///
    /// The degree is stored as non-zero because a finite field always contains
    /// at least the prime field itself.
    pub extension_degree: NonZeroU32,
}

impl FiniteFieldDescriptor {
    /// Builds a descriptor and performs basic consistency checks.
    pub fn new(
        characteristic: impl Into<BigInt>,
        extension_degree: NonZeroU32,
    ) -> Result<Self, FieldError> {
        let characteristic = characteristic.into();
        let Some(characteristic) = characteristic.to_biguint() else {
            return Err(FieldError::InvalidModulus {
                modulus: characteristic.to_string(),
            });
        };
        Self::new_biguint(characteristic, extension_degree)
    }

    /// Builds a descriptor with an arbitrary-precision characteristic.
    pub fn new_biguint(
        characteristic: BigUint,
        extension_degree: NonZeroU32,
    ) -> Result<Self, FieldError> {
        let descriptor = Self {
            characteristic,
            extension_degree,
        };
        descriptor.check()?;
        Ok(descriptor)
    }

    /// Performs lightweight structural checks only.
    pub fn check(&self) -> Result<(), FieldError> {
        if !is_valid_big_field_modulus(&self.characteristic) {
            return Err(FieldError::InvalidModulus {
                modulus: self.characteristic.to_string(),
            });
        }

        Ok(())
    }

    /// Computes the field cardinality as an arbitrary-precision integer.
    pub fn cardinality_biguint(&self) -> BigUint {
        self.characteristic.pow(self.extension_degree.get())
    }

    /// Computes the exact field cardinality.
    pub fn cardinality(&self) -> BigUint {
        self.cardinality_biguint()
    }

    /// Returns a compact educational string such as `F_17` or `F_(43^2)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn pretty(&self) -> String {
        let characteristic = &self.characteristic;
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

    use num_bigint::BigUint;

    use crate::fields::traits::FiniteField;
    use crate::fields::{FieldError, finite_field_descriptor::FiniteFieldDescriptor};

    type F17 = crate::fields::Fp17;

    #[test]
    fn finite_field_descriptor_pretty_formats_prime_fields() {
        let descriptor = FiniteFieldDescriptor::new(17, NonZeroU32::new(1).expect("1 is non-zero"))
            .expect("F17 descriptor should be valid");

        assert_eq!(descriptor.pretty(), "F_17");
        assert_eq!(format!("{descriptor}"), descriptor.pretty());
        assert_eq!(descriptor.cardinality(), BigUint::from(17u8));
    }

    #[test]
    fn finite_field_descriptor_pretty_formats_extension_fields() {
        let descriptor = FiniteFieldDescriptor::new(43, NonZeroU32::new(2).expect("2 is non-zero"))
            .expect("F43^2 descriptor should be valid");

        assert_eq!(descriptor.pretty(), "F_(43^2)");
        assert_eq!(format!("{descriptor}"), descriptor.pretty());
        assert_eq!(descriptor.cardinality(), BigUint::from(43u8).pow(2));
    }

    #[test]
    fn finite_field_descriptor_rejects_invalid_modulus() {
        let error = FiniteFieldDescriptor::new(1, NonZeroU32::new(1).expect("1 is non-zero"))
            .expect_err("characteristic 1 should be rejected");

        assert_eq!(
            error,
            FieldError::InvalidModulus {
                modulus: "1".into()
            }
        );
    }

    #[test]
    fn finite_field_descriptor_keeps_large_cardinality_exact() {
        let characteristic = (BigUint::from(1u8) << 128usize) + BigUint::from(57u8);
        let descriptor = FiniteFieldDescriptor::new_biguint(
            characteristic.clone(),
            NonZeroU32::new(2).expect("2 is non-zero"),
        )
        .expect("large descriptor should be valid");

        assert_eq!(descriptor.cardinality_biguint(), characteristic.pow(2));
    }

    #[test]
    fn finite_field_trait_descriptor_uses_the_metadata_of_the_type() {
        let descriptor =
            F17::descriptor().expect("crate::fields::Fp17 should have a valid descriptor");

        assert_eq!(
            descriptor,
            FiniteFieldDescriptor::new(17, NonZeroU32::new(1).expect("1 is non-zero"))
                .expect("F17 descriptor should be valid")
        );
    }
}
