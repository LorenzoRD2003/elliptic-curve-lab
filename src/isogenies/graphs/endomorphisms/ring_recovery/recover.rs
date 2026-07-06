use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    traits::{CurveModel, FrobeniusTraceCurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, Field, FiniteField, SqrtField};
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphBuilder, IsogenyGraphError, IsogenyGraphNodeId,
    endomorphisms::{
        EndomorphismRingLevelRecoveryError, EndomorphismRingLevelRecoveryReport,
        LocalEndomorphismRingLevelReport,
    },
};

impl<F> IsogenyGraph<ShortWeierstrassCurve<F>>
where
    F: Field + EnumerableFiniteField + SqrtField + FiniteField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
    <ShortWeierstrassCurve<F> as CurveModel>::BaseField: EnumerableFiniteField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + SqrtField<Elem = <ShortWeierstrassCurve<F> as CurveModel>::Elem>
        + FiniteField,
    <ShortWeierstrassCurve<F> as CurveModel>::Point: Clone + Eq + Hash + PartialEq,
    <ShortWeierstrassCurve<F> as GraphCurveModel>::IsomorphismWitness: Clone + fmt::Debug,
    ShortWeierstrassCurve<F>: FrobeniusTraceCurveModel,
{
    /// Recovers the endomorphism ring compatible with local volcano evidence.
    ///
    /// This is the compact user-facing route for Sutherland §3.3. For the
    /// stored representative `E` at `node_id`, Frobenius gives
    /// `Δ_π = v²D_K` and therefore
    /// `ℤ[π] = O_v ⊆ End(E) ≅ O_u ⊆ O_K`, with `u | v`.
    /// For each supplied prime `ℓ`, the method builds the small rooted
    /// `ℓ`-isogeny graph needed to certify the floor distance `δ_ℓ` and then
    /// recovers the local exponent
    ///
    /// `v_ℓ(u) = v_ℓ(v) - δ_ℓ`.
    ///
    /// The returned [`EndomorphismRingLevelRecoveryReport`] is complete only
    /// when the supplied `primes` cover every prime divisor of `v`; otherwise
    /// it exposes the missing primes and leaves the recovered order absent.
    ///
    /// The local floor paths in the returned local reports belong to the
    /// auxiliary rooted `ℓ`-graphs built from the representative curve. The
    /// report's `node_id` still records the node from the original graph whose
    /// curve is being recovered.
    ///
    /// Complexity: one Frobenius/candidate computation for `node_id`; for
    /// each supplied `ℓ`, one small `ℓ`-graph build to depth `v_ℓ(v) + 1`
    /// so the floor candidate is fully expanded, and one floor-distance
    /// search of cost `Θ(δ_ℓ)`.
    pub fn recover_endomorphism_ring_at(
        &self,
        node_id: IsogenyGraphNodeId,
        primes: &[BigUint],
    ) -> Result<EndomorphismRingLevelRecoveryReport, EndomorphismRingLevelRecoveryError> {
        let curve = self
            .node(node_id)
            .ok_or(IsogenyGraphError::MissingSourceNode(node_id))?
            .representative()
            .clone();
        let candidate_set = self.node_endomorphism_candidates(node_id)?;
        let mut local_reports = Vec::with_capacity(primes.len());

        for prime in primes {
            let local_view = candidate_set.local_view_at(prime)?;
            let degree = prime.to_usize().ok_or_else(|| {
                EndomorphismRingLevelRecoveryError::LocalPrimeTooLarge {
                    prime: prime.clone(),
                }
            })?;
            let max_depth = local_view.frobenius_conductor_valuation() as usize + 1;
            let local_graph = IsogenyGraphBuilder::new(curve.clone(), degree)
                .max_depth(max_depth)
                .deduplicate_by_base_field_isomorphism(true)
                .build()?;
            let floor_path = local_graph.find_shortest_floor_path(IsogenyGraphNodeId(0), prime)?;

            local_reports.push(
                LocalEndomorphismRingLevelReport::from_local_view_floor_path_for_node(
                    local_view, node_id, floor_path,
                )?,
            );
        }

        EndomorphismRingLevelRecoveryReport::from_local_reports(candidate_set, local_reports)
    }
}
