use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet,
    quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
};
use crate::numerics::{PositivePrimeError, valuation_biguint};

/// One `ℓ`-local volcanic-level candidate attached to an imaginary quadratic order.
///
/// For an order `O_f = ℤ + f O_K`, its arithmetic `ℓ`-level is
/// `level_ℓ(O_f) = v_ℓ(f)`.
///
/// This value object is intentionally only arithmetic data attached to one
/// candidate order. It does **not** certify that a curve lies at that level in
/// an actual `ℓ`-isogeny volcano.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VolcanoEndomorphismLevelCandidate {
    prime: BigUint,
    order: ImaginaryQuadraticOrder,
    level: u32,
}

impl ImaginaryQuadraticOrder {
    /// Returns the arithmetic `ℓ`-level candidate `level_ℓ(O_f) = v_ℓ(f)`.
    ///
    /// This uses only the conductor of the current order. It is therefore a
    /// local arithmetic invariant of the candidate order, not a certified
    /// statement about the geometric level of a curve in an `ℓ`-volcano.
    ///
    /// Complexity: prime validation is dominated by `num-prime`. After
    /// validation, the implementation performs `Θ(v_ℓ(f))` exact
    /// big-integer divisions.
    pub fn volcanic_level_at(
        &self,
        prime: &BigUint,
    ) -> Result<VolcanoEndomorphismLevelCandidate, PositivePrimeError> {
        VolcanoEndomorphismLevelCandidate::new(
            prime.clone(),
            self.clone(),
            valuation_biguint(self.conductor(), prime)?,
        )
    }
}

impl EndomorphismRingCandidateSet {
    /// Returns the arithmetic `ℓ`-level candidates for every Frobenius-compatible order `O_f`.
    ///
    /// If the candidate set stores all orders with `ℤ[π] ⊆ O_f ⊆ O_K`, this
    /// annotates each one with its local arithmetic level `v_ℓ(f)`.
    ///
    /// Complexity: dominated by `num-prime`, plus `Θ(τ(v) · v_ℓ(v))`
    /// exact big-integer divisions in the current straightforward pass.
    pub fn volcanic_level_candidates_at(
        &self,
        prime: &BigUint,
    ) -> Result<Vec<VolcanoEndomorphismLevelCandidate>, PositivePrimeError> {
        self.candidate_orders()
            .iter()
            .map(|order| order.volcanic_level_at(prime))
            .collect()
    }
}

impl VolcanoEndomorphismLevelCandidate {
    fn new(
        prime: BigUint,
        order: ImaginaryQuadraticOrder,
        level: u32,
    ) -> Result<Self, PositivePrimeError> {
        valuation_biguint(&prime, &prime)?;

        Ok(Self {
            prime,
            order,
            level,
        })
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the candidate order `O_f`.
    pub fn order(&self) -> &ImaginaryQuadraticOrder {
        &self.order
    }

    /// Returns the arithmetic local level `v_ℓ(f)`.
    pub fn level(&self) -> u32 {
        self.level
    }

    /// Returns the conductor `f` of the candidate order.
    pub fn conductor(&self) -> &BigUint {
        self.order.conductor()
    }

    /// Returns the discriminant `Δ(O_f) = f^2 D_K` of the candidate order.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        self.order.discriminant()
    }
}
