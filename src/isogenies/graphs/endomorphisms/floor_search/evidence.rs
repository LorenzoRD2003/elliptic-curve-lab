use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::One;

use crate::fields::traits::Field;
use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId, endomorphisms::VolcanoSearchError,
};
use crate::numerics::validate_positive_prime;

/// Local status inferred from Proposition 10 for an ordinary `ℓ`-volcano.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum VolcanoFloorStatus {
    /// The expanded node has `deg(v) ≤ 2`, so the ordinary-volcano criterion
    /// places it on the floor.
    OnFloor,
    /// The expanded node has `deg(v) = ℓ + 1`, so the ordinary-volcano
    /// criterion says it is not on the floor.
    NotOnFloor,
    /// The node is present but was not fully expanded, so the observed degree
    /// is only a lower bound.
    UnknownBecausePartialGraph,
    /// The node has `j = 0` or `j = 1728`, where the clean volcano statement
    /// needs special automorphism corrections.
    SpecialJInvariant,
    /// The complete observed degree is neither floor-like nor `ℓ + 1`.
    InconsistentWithVolcanoModel,
}

/// Local floor evidence for one vertex in a stored `ℓ`-isogeny graph.
///
/// The key mathematical convention is Sutherland's Proposition 10:
/// in an ordinary `ℓ`-volcano, a vertex with `deg(v) ≤ 2` is on the floor,
/// while a non-floor vertex has `deg(v) = ℓ + 1`.
///
/// This report uses the stored outgoing edge count as the observed degree,
/// with multiplicity. The status is decisive only when the graph says the node
/// was fully expanded; otherwise it records
/// [`VolcanoFloorStatus::UnknownBecausePartialGraph`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct VolcanoFloorEvidence {
    observed_out_degree: usize,
    status: VolcanoFloorStatus,
}

impl VolcanoFloorEvidence {
    pub(crate) fn new(observed_out_degree: usize, status: VolcanoFloorStatus) -> Self {
        Self {
            observed_out_degree,
            status,
        }
    }

    /// Returns the stored outgoing degree, counted with edge multiplicity.
    pub(crate) fn observed_out_degree(&self) -> usize {
        self.observed_out_degree
    }

    /// Returns the local floor status.
    pub(crate) fn status(&self) -> VolcanoFloorStatus {
        self.status
    }
}

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Reads the local floor evidence at one graph node.
    ///
    /// This is the smallest public surface for §3.1-style floor detection. It
    /// does not compute endomorphism rings and does not prove that the ambient
    /// component is ordinary. It only applies the local degree test when the
    /// graph builder has fully expanded the requested node.
    pub(crate) fn floor_evidence_at(
        &self,
        node_id: IsogenyGraphNodeId,
        prime: &BigUint,
    ) -> Result<VolcanoFloorEvidence, VolcanoSearchError> {
        validate_positive_prime(prime)?;
        let node = self
            .node(node_id)
            .ok_or(VolcanoSearchError::NodeNotFound { node_id })?;

        let expected_non_floor_degree = prime + BigUint::one();
        let observed_out_degree = self.out_degree(node_id);
        let status = if is_special_j_invariant::<C>(&node.j_invariant()) {
            VolcanoFloorStatus::SpecialJInvariant
        } else if self.node_is_fully_expanded(node_id) != Some(true) {
            VolcanoFloorStatus::UnknownBecausePartialGraph
        } else if observed_out_degree <= 2 {
            VolcanoFloorStatus::OnFloor
        } else if BigUint::from(observed_out_degree) == expected_non_floor_degree {
            VolcanoFloorStatus::NotOnFloor
        } else {
            VolcanoFloorStatus::InconsistentWithVolcanoModel
        };

        Ok(VolcanoFloorEvidence::new(observed_out_degree, status))
    }
}

fn is_special_j_invariant<C: GraphCurveModel>(j_invariant: &C::Elem) -> bool {
    C::BaseField::eq(j_invariant, &C::BaseField::zero())
        || C::BaseField::eq(j_invariant, &C::BaseField::from_i64(1728))
}
