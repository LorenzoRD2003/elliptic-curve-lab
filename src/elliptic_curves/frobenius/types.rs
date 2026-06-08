use crate::fields::{FiniteField, FiniteFieldDescriptor};

/// Metadata for the absolute Frobenius `π_p^k`.
///
/// This value object records only the characteristic `p` and the iterate `k`.
/// It does not yet carry a separate geometric codomain curve, so callers
/// should read it as descriptive metadata rather than as a full curve-map
/// witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AbsoluteFrobenius {
    /// Prime characteristic `p`.
    pub characteristic: u64,
    /// Iterate count in `π_p^power`.
    pub power: u32,
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
    ///
    /// This helper records the prime characteristic only; it does not attempt
    /// to decide whether the resulting map is already an endomorphism of a
    /// chosen curve model.
    pub fn for_field<F: FiniteField>(power: u32) -> Self {
        Self::new(F::characteristic(), power)
    }

    /// Returns whether this Frobenius iterate is the identity map.
    pub fn is_identity(&self) -> bool {
        self.power == 0
    }
}

/// Metadata for the relative Frobenius `π_q^k` on a finite base field.
///
/// The stored field descriptor makes the choice of `q = p^r` explicit, which
/// is useful later when contrasting:
///
/// - absolute Frobenius `π_p` on geometric points
/// - relative Frobenius `π_q` as an endomorphism of a curve over `F_q`
///
/// In particular, if a point is fixed by this map for `power = 1`, it is a
/// natural candidate for being rational over the represented base field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelativeFrobenius {
    /// Finite base field metadata for `F_q`.
    pub base_field: FiniteFieldDescriptor,
    /// Iterate count in `π_q^power`.
    pub power: u32,
}

impl RelativeFrobenius {
    /// Builds relative Frobenius metadata from explicit parameters.
    pub fn new(base_field: FiniteFieldDescriptor, power: u32) -> Self {
        Self { base_field, power }
    }

    /// Builds relative Frobenius metadata for a concrete finite field family.
    ///
    /// If `F` has size `q = p^r`, the returned metadata represents `π_q^power`.
    pub fn for_field<F: FiniteField>(power: u32) -> Self {
        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .expect("finite field implementations should expose internally consistent metadata");

        Self::new(base_field, power)
    }

    /// Returns whether this Frobenius iterate is the identity map.
    pub fn is_identity(&self) -> bool {
        self.power == 0
    }
}
