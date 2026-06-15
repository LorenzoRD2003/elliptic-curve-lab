use num_bigint::BigUint;
use std::ops::RangeInclusive;

use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet,
    quadratic_orders::QuadraticDiscriminantFactorization,
};
use crate::numerics::{PositivePrimeError, valuation_biguint};

/// Local `ℓ`-adic view of the conductor gap between `ℤ[π]` and `O_K`.
///
/// If the Frobenius discriminant factors as `Δ_π = v^2 D_K`, then `ℤ[π] = O_v`
/// and the `ℓ`-local ambiguity is measured by`v_ℓ(v)`. This value object
/// records exactly that exponent for one chosen prime `ℓ`.
///
/// The candidate local conductors then have exponents `0 <= b <= v_ℓ(v)`,
/// interpolating between the maximal order `O_K = O_1` at `b = 0` and the
/// Frobenius order `ℤ[π] = O_v` at `b = v_ℓ(v)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndomorphismRingLocalView {
    prime: BigUint,
    frobenius_conductor_valuation: u32,
}

impl QuadraticDiscriminantFactorization {
    /// Returns the local `ℓ`-adic view determined by `v_ℓ(v)`.
    ///
    /// If `Δ_π = v^2 D_K`, this measures the `ℓ`-primary part of the
    /// conductor of `ℤ[π] = O_v`.
    ///
    /// Complexity: prime validation is dominated by `num-prime`. After validation,
    /// the implementation performs `Θ(v_ℓ(v))` exact big-integer divisions.
    pub fn local_view_at(
        &self,
        prime: &BigUint,
    ) -> Result<EndomorphismRingLocalView, PositivePrimeError> {
        EndomorphismRingLocalView::new(prime.clone(), valuation_biguint(self.conductor(), prime)?)
    }
}

impl EndomorphismRingCandidateSet {
    /// Returns the local `ℓ`-adic view determined by the Frobenius conductor `v`.
    ///
    /// Complexity: prime validation is dominated by `num-prime`. After validation,
    /// the implementation performs `Θ(v_ℓ(v))` exact big-integer divisions.
    pub fn local_view_at(
        &self,
        prime: &BigUint,
    ) -> Result<EndomorphismRingLocalView, PositivePrimeError> {
        self.factorization().local_view_at(prime)
    }
}

impl EndomorphismRingLocalView {
    fn new(prime: BigUint, frobenius_conductor_valuation: u32) -> Result<Self, PositivePrimeError> {
        valuation_biguint(&prime, &prime)?;

        Ok(Self {
            prime,
            frobenius_conductor_valuation,
        })
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the exponent `v_ℓ(v)` of the Frobenius conductor.
    pub fn frobenius_conductor_valuation(&self) -> u32 {
        self.frobenius_conductor_valuation
    }

    /// Returns whether the `ℓ`-local conductor gap is trivial, i.e.,
    /// whether `ℓ` does not divide the conductor of `ℤ[π]`.
    pub fn is_trivial(&self) -> bool {
        self.frobenius_conductor_valuation == 0
    }

    /// Returns the possible local conductor exponents `0 <= b <= v_ℓ(v)`.
    ///
    /// This is the local chain of candidate exponents between `O_K` and `ℤ[π]`.
    pub fn local_conductor_exponents(&self) -> RangeInclusive<u32> {
        0..=self.frobenius_conductor_valuation
    }
}
