use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::One;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId,
    endomorphisms::{VolcanoSearchError, floor_search::evidence::VolcanoFloorStatus},
};
use crate::numerics::validate_positive_prime;

/// Report produced by the `FindShortestPathToFloor`-style search.
///
/// For a complete ordinary `ℓ`-volcano, Sutherland §3.1 uses the fact that
/// from any non-floor vertex at most two adjacent vertices are not descending.
/// Starting three non-backtracking paths therefore guarantees that at least one
/// path descends toward the floor. This report records the first floor path
/// certified by that parallel search, so [`Self::distance_to_floor`] is the
/// volcano distance
///
/// `δ(v) = dist(v, V_d)`.
///
/// Algorithmically, the search keeps at most three active paths, so its
/// abstract running time is `Θ(δ(v))`, independent of `ℓ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShortestFloorPathReport {
    prime: BigUint,
    start: IsogenyGraphNodeId,
    floor: IsogenyGraphNodeId,
    path: Vec<IsogenyGraphNodeId>,
}

impl ShortestFloorPathReport {
    pub(crate) fn new(
        prime: BigUint,
        start: IsogenyGraphNodeId,
        floor: IsogenyGraphNodeId,
        path: Vec<IsogenyGraphNodeId>,
    ) -> Self {
        Self {
            prime,
            start,
            floor,
            path,
        }
    }

    /// Returns the chosen prime `ℓ`.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the start vertex.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.start
    }

    /// Returns the certified floor vertex.
    pub fn floor(&self) -> IsogenyGraphNodeId {
        self.floor
    }

    /// Returns the certified shortest path from `start` to the floor.
    pub fn path(&self) -> &[IsogenyGraphNodeId] {
        &self.path
    }

    /// Returns the volcano distance `δ(v) = dist(v, V_d)`.
    pub fn distance_to_floor(&self) -> usize {
        self.path.len().saturating_sub(1)
    }

    /// Computes the volcano level when the total depth `d` is already known.
    ///
    /// The convention is `level(v) = d - δ(v)`, where level `0` is the
    /// surface and level `d` is the floor. Returns `None` if this report's
    /// distance is larger than the supplied total depth.
    pub fn level_from_total_depth(&self, depth: usize) -> Option<usize> {
        depth.checked_sub(self.distance_to_floor())
    }
}

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Finds a shortest path from `start` to the floor of an ordinary
    /// `ℓ`-volcano.
    ///
    /// This is the graph-side `FindShortestPathToFloor` milestone from
    /// Sutherland §3.1. If `start` is not already on the floor, the algorithm
    /// starts up to three non-backtracking paths from distinct outgoing
    /// neighbors and advances them in parallel. In a complete ordinary
    /// `ℓ`-volcano, at most two adjacent vertices are non-descending, so at
    /// least one of these paths descends monotonically to the floor.
    ///
    /// The method is deliberately conservative about graph evidence. It only
    /// certifies `δ(v)` from complete local degree data; partial boundary
    /// nodes, special `j`-invariants, or degrees inconsistent with the volcano
    /// model are reported as errors rather than converted into a misleading
    /// shortest-path claim.
    ///
    /// Complexity: `Θ(δ(v))`, because the volcano argument keeps the search
    /// width bounded by three paths, independent of `ℓ`.
    pub fn find_shortest_floor_path(
        &self,
        start: IsogenyGraphNodeId,
        prime: &BigUint,
    ) -> Result<ShortestFloorPathReport, VolcanoSearchError> {
        validate_positive_prime(prime)?;

        let start_evidence = self.floor_evidence_at(start, prime)?;
        match start_evidence.status() {
            VolcanoFloorStatus::OnFloor => {
                return Ok(ShortestFloorPathReport::new(
                    prime.clone(),
                    start,
                    start,
                    vec![start],
                ));
            }
            VolcanoFloorStatus::NotOnFloor => {}
            VolcanoFloorStatus::UnknownBecausePartialGraph => {
                return Err(VolcanoSearchError::NodeNotFullyExpanded { node_id: start });
            }
            VolcanoFloorStatus::SpecialJInvariant => {
                return Err(VolcanoSearchError::SpecialJInvariant { node_id: start });
            }
            VolcanoFloorStatus::InconsistentWithVolcanoModel => {
                return Err(VolcanoSearchError::InconsistentWithVolcanoModel {
                    node_id: start,
                    observed_out_degree: start_evidence.observed_out_degree(),
                    expected_non_floor_degree: prime + BigUint::one(),
                });
            }
        }

        let mut paths = self
            .non_backtracking_outgoing_neighbors(start, None)
            .into_iter()
            .take(3)
            .map(|neighbor| TrackedFloorPath::from_first_step(start, neighbor))
            .collect::<Vec<_>>();

        if paths.is_empty() {
            return Err(VolcanoSearchError::NoNonBacktrackingNeighbor { node_id: start });
        }

        loop {
            for path in &paths {
                let current = path.current();
                let evidence = self.floor_evidence_at(current, prime)?;
                match evidence.status() {
                    VolcanoFloorStatus::OnFloor => {
                        return Ok(ShortestFloorPathReport::new(
                            prime.clone(),
                            start,
                            current,
                            path.nodes.clone(),
                        ));
                    }
                    VolcanoFloorStatus::NotOnFloor => {}
                    VolcanoFloorStatus::UnknownBecausePartialGraph => {
                        return Err(VolcanoSearchError::NodeNotFullyExpanded { node_id: current });
                    }
                    VolcanoFloorStatus::SpecialJInvariant => {
                        return Err(VolcanoSearchError::SpecialJInvariant { node_id: current });
                    }
                    VolcanoFloorStatus::InconsistentWithVolcanoModel => {
                        return Err(VolcanoSearchError::InconsistentWithVolcanoModel {
                            node_id: current,
                            observed_out_degree: evidence.observed_out_degree(),
                            expected_non_floor_degree: prime + BigUint::one(),
                        });
                    }
                }
            }

            let mut extended_paths = Vec::with_capacity(paths.len());
            for mut path in paths {
                let Some(next) = self
                    .non_backtracking_outgoing_neighbors(path.current(), Some(path.previous()))
                    .into_iter()
                    .next()
                else {
                    continue;
                };

                if path.try_push(next) {
                    extended_paths.push(path);
                }
            }

            if extended_paths.is_empty() {
                return Err(VolcanoSearchError::NoFloorPathFound { start });
            }

            paths = extended_paths;
        }
    }
}

#[derive(Clone, Debug)]
struct TrackedFloorPath {
    nodes: Vec<IsogenyGraphNodeId>,
    visited: HashSet<IsogenyGraphNodeId>,
}

impl TrackedFloorPath {
    fn from_first_step(start: IsogenyGraphNodeId, next: IsogenyGraphNodeId) -> Self {
        Self {
            nodes: vec![start, next],
            visited: HashSet::from([start, next]),
        }
    }

    fn current(&self) -> IsogenyGraphNodeId {
        *self
            .nodes
            .last()
            .expect("tracked floor paths are never empty")
    }

    fn previous(&self) -> IsogenyGraphNodeId {
        self.nodes[self.nodes.len() - 2]
    }

    fn try_push(&mut self, next: IsogenyGraphNodeId) -> bool {
        if !self.visited.insert(next) {
            return false;
        }

        self.nodes.push(next);
        true
    }
}
