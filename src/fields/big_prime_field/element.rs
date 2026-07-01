use core::fmt;

use num_bigint::BigUint;

/// Canonical residue class in one runtime prime field.
///
/// Relative to an ambient field `F_p`, values are stored in the canonical
/// range `0 ≤ x < p`.
///
/// The element itself intentionally does not carry a runtime modulus. The
/// ambient [`super::BigPrimeField`] owns that context and is responsible for
/// arithmetic, equality, and inversion.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BigPrimeFieldElem(BigUint);

impl BigPrimeFieldElem {
    /// Builds one canonical residue class from a value already reduced modulo
    /// the ambient prime.
    pub(crate) fn new_canonical(value: BigUint) -> Self {
        Self(value)
    }

    /// Returns the stored canonical representative in `ℤ`.
    pub fn value(&self) -> &BigUint {
        &self.0
    }
}

impl fmt::Display for BigPrimeFieldElem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
