use core::fmt;

use num_bigint::BigUint;

/// Characteristic of a field family.
///
/// A field has either characteristic `0` or positive prime characteristic `p`.
/// Positive characteristic is stored as a `BigUint` so future static large
/// prime backends can keep the base trait honest instead of truncating `p` to
/// `u64`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FieldCharacteristic {
    /// Characteristic zero.
    Zero,
    /// Positive prime characteristic `p`.
    Positive(BigUint),
}

impl FieldCharacteristic {
    /// Returns whether the characteristic is zero.
    pub fn is_zero(&self) -> bool {
        matches!(self, Self::Zero)
    }

    /// Returns the positive characteristic as a `BigUint`, or `None` in
    /// characteristic zero.
    pub fn to_positive_biguint(&self) -> Option<BigUint> {
        match self {
            Self::Zero => None,
            Self::Positive(characteristic) => Some(characteristic.clone()),
        }
    }

    /// Returns the characteristic as an exact non-negative integer.
    pub fn to_biguint(&self) -> BigUint {
        match self {
            Self::Zero => BigUint::from(0u8),
            Self::Positive(characteristic) => characteristic.clone(),
        }
    }
}

impl fmt::Display for FieldCharacteristic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => f.write_str("0"),
            Self::Positive(characteristic) => write!(f, "{characteristic}"),
        }
    }
}
