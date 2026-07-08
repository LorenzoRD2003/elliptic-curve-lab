use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{
        PrimeNormIdealError, ramified_prime_ideal::RamifiedPrimeIdeal,
        split_prime_ideal::SplitPrimeIdeal,
    },
    quadratic_orders::ImaginaryQuadraticOrder,
};

/// A supported prime-norm ideal in an imaginary quadratic order.
///
/// This public type is intentionally opaque for now. It gives downstream code
/// a stable place to accept “a supported ideal of prime norm” without exposing
/// the concrete split or ramified representations, or implying that inert
/// primes, conductor-dividing primes, multiplication, or class-group operations
/// already exist.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrimeNormIdeal {
    kind: PrimeNormIdealKind,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PrimeNormIdealKind {
    Split(SplitPrimeIdeal),
    Ramified(RamifiedPrimeIdeal),
}

impl PrimeNormIdeal {
    /// Builds a split prime ideal of norm `ℓ`.
    ///
    /// This is a convenience wrapper around [`SplitPrimeIdeal::new`].
    pub fn split(
        order: ImaginaryQuadraticOrder,
        ell: BigUint,
        root: BigUint,
    ) -> Result<Self, PrimeNormIdealError> {
        Ok(Self {
            kind: PrimeNormIdealKind::Split(SplitPrimeIdeal::new(order, ell, root)?),
        })
    }

    /// Builds the unique ramified prime ideal of norm `ℓ`.
    ///
    /// The repeated root modulo `ℓ` is derived from
    /// [`ImaginaryQuadraticOrder::prime_behavior`]; callers do not choose it.
    pub fn ramified(
        order: ImaginaryQuadraticOrder,
        ell: BigUint,
    ) -> Result<Self, PrimeNormIdealError> {
        Ok(Self {
            kind: PrimeNormIdealKind::Ramified(RamifiedPrimeIdeal::new(order, ell)?),
        })
    }

    /// Returns the imaginary quadratic order containing the ideal.
    pub fn order(&self) -> &ImaginaryQuadraticOrder {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.order(),
            PrimeNormIdealKind::Ramified(ideal) => ideal.order(),
        }
    }

    /// Returns the prime norm `ℓ`.
    pub fn norm(&self) -> &BigUint {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.norm(),
            PrimeNormIdealKind::Ramified(ideal) => ideal.norm(),
        }
    }

    /// Returns the local root of `Δ` modulo `ℓ` carried by this ideal.
    ///
    /// In the split case this is the selected root; in the ramified case this
    /// is the repeated root.
    pub fn root_mod_ell(&self) -> &BigUint {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.root(),
            PrimeNormIdealKind::Ramified(ideal) => ideal.root(),
        }
    }

    /// Returns whether this ideal is represented by the split-prime family.
    pub fn is_split(&self) -> bool {
        matches!(self.kind, PrimeNormIdealKind::Split(_))
    }

    /// Returns whether this ideal is represented by the ramified-prime family.
    pub fn is_ramified(&self) -> bool {
        matches!(self.kind, PrimeNormIdealKind::Ramified(_))
    }

    /// Returns the conjugate prime-norm ideal.
    pub fn conjugate(&self) -> Self {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => Self {
                kind: PrimeNormIdealKind::Split(ideal.conjugate()),
            },
            PrimeNormIdealKind::Ramified(ideal) => Self {
                kind: PrimeNormIdealKind::Ramified(ideal.conjugate()),
            },
        }
    }
}
