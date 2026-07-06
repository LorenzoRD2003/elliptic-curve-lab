use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::elliptic_curves::{
    endomorphisms::candidate_sets::EndomorphismRingLocalView, traits::FrobeniusTraceCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId,
    endomorphisms::{EndomorphismRingLevelRecoveryError, ShortestFloorPathReport},
};

/// Local recovery report for the conductor exponent of `End(E)` at one prime `ℓ`.
///
/// For an ordinary elliptic curve over `F_q`, Frobenius satisfies
/// `π² - tπ + q = 0`, and the Frobenius discriminant is
/// `Δ_π = t² - 4q`. If `Δ_π = v²D_K`, then `ℤ[π] = O_v` and
/// `[O_K : ℤ[π]] = v`. Since
///
/// `ℤ[π] ⊆ End(E) ⊆ O_K`,
///
/// there is a divisor `u | v` such that `End(E) ≅ O_u`. For a prime-power
/// factor `ℓ^e || v`, Sutherland §3.3 relates the distance to the floor in the
/// `ℓ`-volcano to the local conductor exponent:
///
/// `δ = e - d`, where `d = v_ℓ(u)`.
///
/// Equivalently, this report recovers `v_ℓ(u) = e - δ`.
///
/// The report is *local*: it identifies only the `ℓ`-part of the
/// endomorphism-ring conductor. Recovering the full conductor requires the
/// corresponding data for every prime divisor of `v`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalEndomorphismRingLevelReport {
    local_view: EndomorphismRingLocalView,
    floor_path: ShortestFloorPathReport,
}

impl LocalEndomorphismRingLevelReport {
    pub(crate) fn from_local_view_and_floor_path(
        local_view: EndomorphismRingLocalView,
        floor_path: ShortestFloorPathReport,
    ) -> Result<Self, EndomorphismRingLevelRecoveryError> {
        let distance_to_floor = floor_path.distance_to_floor();
        if distance_to_floor > local_view.frobenius_conductor_valuation() as usize {
            return Err(
                EndomorphismRingLevelRecoveryError::DistanceExceedsFrobeniusConductorValuation {
                    node_id: floor_path.start(),
                    prime: local_view.prime().clone(),
                    distance_to_floor,
                    frobenius_conductor_valuation: local_view.frobenius_conductor_valuation(),
                },
            );
        }

        Ok(Self {
            local_view,
            floor_path,
        })
    }

    /// Returns the node whose local ring level was recovered.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.floor_path.start()
    }

    /// Returns the chosen local prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        self.local_view.prime()
    }

    /// Returns `e = v_ℓ(v)`, where `Δ_π = v²D_K`.
    pub fn frobenius_conductor_valuation(&self) -> u32 {
        self.local_view.frobenius_conductor_valuation()
    }

    /// Returns the certified volcano distance `δ = dist(E, V_d)`.
    pub fn distance_to_floor(&self) -> usize {
        self.floor_path.distance_to_floor()
    }

    /// Returns the recovered local conductor exponent `d = v_ℓ(u)`.
    ///
    /// If `End(E) ≅ O_u` and `ℓᵉ || v`, then Sutherland's volcano relation
    /// gives `d = e - δ`, where `δ` is the distance from the node to the floor.
    pub fn recovered_conductor_valuation(&self) -> u32 {
        self.frobenius_conductor_valuation() - self.distance_to_floor() as u32
    }

    /// Returns the shortest floor path certifying `δ`.
    pub fn floor_path(&self) -> &ShortestFloorPathReport {
        &self.floor_path
    }
}

impl<C: GraphCurveModel + FrobeniusTraceCurveModel> IsogenyGraph<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem> + FiniteField,
    C::Point: Clone + Eq + Hash + PartialEq,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Recovers the local conductor exponent `v_ℓ(u)` for `End(E) ≅ O_u`.
    ///
    /// This is the local ring-recovery step from Sutherland §3.3. The method:
    ///
    /// 1. derives the Frobenius-compatible candidate orders from the node,
    /// 2. reads `e = v_ℓ(v)` from `Δ_π = v²D_K`,
    /// 3. certifies the shortest distance `δ` from the node to the floor of
    ///    the observed `ℓ`-volcano,
    /// 4. returns the local exponent `d = e - δ`.
    ///
    /// The result is only one local component of the conductor of `End(E)`.
    /// It does not by itself recover the full order `O_u` unless every prime
    /// divisor of `v` has also been recovered.
    ///
    /// Complexity: one exhaustive Frobenius-trace computation for the stored
    /// curve, plus `num-prime` dominated conductor arithmetic, plus `Θ(δ)` for
    /// the shortest floor search.
    pub fn recover_endomorphism_ring_level_at(
        &self,
        node_id: IsogenyGraphNodeId,
        prime: &BigUint,
    ) -> Result<LocalEndomorphismRingLevelReport, EndomorphismRingLevelRecoveryError> {
        let candidate_set = self.node_endomorphism_candidates(node_id)?;
        let local_view = candidate_set.local_view_at(prime)?;
        let floor_path = self.find_shortest_floor_path(node_id, prime)?;

        LocalEndomorphismRingLevelReport::from_local_view_and_floor_path(local_view, floor_path)
    }
}
