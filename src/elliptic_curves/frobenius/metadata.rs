use crate::fields::{finite_field_descriptor::FiniteFieldDescriptor, traits::FiniteField};

/// Metadata for the absolute Frobenius `π_p^k`.
///
/// This value object records only the characteristic `p` and the iterate `k`.
/// It does not yet carry a separate geometric codomain curve, so callers
/// should read it as descriptive metadata rather than as a full curve-map
/// witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AbsoluteFrobenius {
    characteristic: u64,
    power: u32,
}

impl AbsoluteFrobenius {
    /// Builds absolute Frobenius metadata from explicit parameters.
    pub fn new(characteristic: u64, power: u32) -> Self {
        Self {
            characteristic,
            power,
        }
    }

    /// Builds absolute Frobenius metadata for a concrete finite field family.
    pub fn for_field<F: FiniteField>(power: u32) -> Self {
        Self::new(F::characteristic(), power)
    }

    /// Returns whether this Frobenius iterate is the identity map.
    pub fn is_identity(&self) -> bool {
        self.power == 0
    }

    pub fn characteristic(&self) -> u64 {
        self.characteristic
    }

    pub fn power(&self) -> u32 {
        self.power
    }
}

/// Metadata for the relative Frobenius `π_q^k` on a finite base field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelativeFrobenius {
    base_field: FiniteFieldDescriptor,
    power: u32,
}

impl RelativeFrobenius {
    /// Builds relative Frobenius metadata from explicit parameters.
    pub fn new(base_field: FiniteFieldDescriptor, power: u32) -> Self {
        Self { base_field, power }
    }

    /// Builds relative Frobenius metadata for a concrete finite field family.
    pub fn for_field<F: FiniteField>(power: u32) -> Self {
        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .expect("finite field implementations should expose internally consistent metadata");

        Self::new(base_field, power)
    }

    /// Returns whether this Frobenius iterate is the identity map.
    pub fn is_identity(&self) -> bool {
        self.power == 0
    }

    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    pub fn power(&self) -> u32 {
        self.power
    }
}
