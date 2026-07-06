use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId,
    endomorphisms::IsogenyEdgeEndomorphismTentativeRelation,
};

/// Cached graph-side altimeter evidence for an ordinary `ℓ`-volcano.
///
/// The altimeter records the certified distance
///
/// `δ(v) = dist(v, V_d)`
///
/// from each node to the floor when [`IsogenyGraph::find_shortest_floor_path`]
/// can certify it from complete local degree evidence. Missing entries mean
/// the graph did not provide enough ordinary-volcano evidence for that node.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct VolcanoAltimeterEvidence {
    distances_to_floor: Vec<Option<usize>>,
}

impl VolcanoAltimeterEvidence {
    /// Builds cached altimeter evidence for every stored node.
    ///
    /// A failed shortest-floor search is intentionally degraded to `None`.
    /// The altimeter is optional graph evidence for refinement, not a
    /// diagnostic API; callers should fall back to weaker evidence when a
    /// node's `δ(v)` cannot be certified.
    ///
    /// Complexity: `Θ(Σ_v δ(v))`, equivalently `O(n δ_max)` for `n` stored
    /// nodes and maximum certified floor distance `δ_max`.
    pub(crate) fn from_graph<C: GraphCurveModel>(graph: &IsogenyGraph<C>, prime: &BigUint) -> Self
    where
        C::Point: Clone + Eq + Hash,
        C::IsomorphismWitness: Clone + fmt::Debug,
    {
        let distances_to_floor = graph
            .nodes()
            .iter()
            .map(|node| {
                graph
                    .find_shortest_floor_path(node.id(), prime)
                    .ok()
                    .map(|report| report.distance_to_floor())
            })
            .collect();

        Self { distances_to_floor }
    }

    pub(crate) fn distance_to_floor(&self, node_id: IsogenyGraphNodeId) -> Option<usize> {
        self.distances_to_floor.get(node_id.0).copied().flatten()
    }

    /// Classifies one directed edge by comparing endpoint altitudes.
    ///
    /// If both endpoint distances are certified, then:
    ///
    /// - equal distances mean horizontal;
    /// - source distance one larger than target distance means descending;
    /// - target distance one larger than source distance means ascending.
    ///
    /// Jumps larger than one level are rejected as incompatible with the local
    /// ordinary-volcano edge shape represented by this evidence.
    pub(crate) fn relation_for(
        &self,
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
    ) -> Option<IsogenyEdgeEndomorphismTentativeRelation> {
        let source_distance = self.distance_to_floor(source)?;
        let target_distance = self.distance_to_floor(target)?;

        IsogenyEdgeEndomorphismTentativeRelation::from_floor_distances(
            source_distance,
            target_distance,
        )
    }
}

#[cfg(test)]
mod tests;
