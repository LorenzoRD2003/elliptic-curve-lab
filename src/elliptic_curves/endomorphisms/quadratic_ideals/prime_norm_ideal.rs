use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{PrimeNormIdealError, split_prime_ideal::SplitPrimeIdeal},
    quadratic_orders::ImaginaryQuadraticOrder,
};

/// A supported prime-norm ideal in an imaginary quadratic order.
///
/// This public type is intentionally opaque for now. It gives downstream code
/// a stable place to accept “a supported ideal of prime norm” without exposing
/// the current one-variant representation or implying that inert primes,
/// conductor-dividing primes, ramified ideals, multiplication, or class-group
/// operations already exist.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrimeNormIdeal {
    kind: PrimeNormIdealKind,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PrimeNormIdealKind {
    Split(SplitPrimeIdeal),
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

    /// Returns the imaginary quadratic order containing the ideal.
    pub fn order(&self) -> &ImaginaryQuadraticOrder {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.order(),
        }
    }

    /// Returns the prime norm `ℓ`.
    pub fn norm(&self) -> &BigUint {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.norm(),
        }
    }

    /// Returns the selected split root when this ideal is represented by a
    /// split prime above `ℓ`.
    ///
    /// This accessor is intentionally observational: callers may inspect the
    /// current split-root witness, but the concrete split-ideal type remains
    /// crate-internal.
    pub fn split_root(&self) -> &BigUint {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => ideal.root(),
        }
    }

    /// Returns the conjugate prime-norm ideal.
    pub fn conjugate(&self) -> Self {
        match &self.kind {
            PrimeNormIdealKind::Split(ideal) => Self {
                kind: PrimeNormIdealKind::Split(ideal.conjugate()),
            },
        }
    }
}
