use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::One;

use crate::isogenies::graphs::{
    GraphCurveModel, IsogenyGraph, IsogenyGraphNodeId,
    endomorphisms::floor_search::VolcanoSearchError,
};
use crate::numerics::validate_positive_prime;

use super::evidence::VolcanoFloorStatus;

/// Report produced by the first `FindFloor`-style walk.
///
/// This report records the path actually followed by the selected strategy.
/// It is not a shortest-path report and therefore does not yet certify the
/// minimal distance `δ(E)` from the start vertex to the floor. That stronger
/// surface belongs to a later `FindShortestPathToFloor` milestone.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FloorPathReport {
    prime: BigUint,
    start: IsogenyGraphNodeId,
    floor: IsogenyGraphNodeId,
    path: Vec<IsogenyGraphNodeId>,
}

impl FloorPathReport {
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

    /// Returns the first vertex on the followed path with floor evidence.
    pub fn floor(&self) -> IsogenyGraphNodeId {
        self.floor
    }

    /// Returns the followed non-backtracking path.
    pub fn path(&self) -> &[IsogenyGraphNodeId] {
        &self.path
    }

    /// Returns the length of the followed path.
    ///
    /// This is a `FindFloor` path length, not yet the minimal distance `δ(E)`.
    pub fn distance_to_floor(&self) -> usize {
        self.path.len().saturating_sub(1)
    }
}

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Follows a first `FindFloor`-style path until local floor evidence appears.
    ///
    /// This method implements the deliberately small milestone based on
    /// Sutherland §3.1. It uses local complete outgoing degree evidence:
    ///
    /// - `deg(v) ≤ 2` means the current vertex is on the floor;
    /// - `deg(v) = ℓ + 1` means the current vertex is not on the floor.
    ///
    /// The current walk is deterministic and reproducible: it follows the
    /// first stored outgoing edge that does not immediately backtrack. It is
    /// not the randomized algorithm from the paper and it is not
    /// `FindShortestPathToFloor`, so [`FloorPathReport::distance_to_floor`]
    /// records only the path length found by this deterministic walk.
    pub fn find_floor_path(
        &self,
        start: IsogenyGraphNodeId,
        prime: &BigUint,
    ) -> Result<FloorPathReport, VolcanoSearchError> {
        validate_positive_prime(prime)?;

        let mut path = vec![start];
        let mut visited = HashSet::from([start]);
        let mut previous = None;
        let mut current = start;

        loop {
            let evidence = self.floor_evidence_at(current, prime)?;
            match evidence.status() {
                VolcanoFloorStatus::OnFloor => {
                    return Ok(FloorPathReport::new(prime.clone(), start, current, path));
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

            let next = self
                .first_non_backtracking_outgoing_neighbor(current, previous)
                .ok_or(VolcanoSearchError::NoNonBacktrackingNeighbor { node_id: current })?;

            if !visited.insert(next) {
                return Err(VolcanoSearchError::CycleDetectedBeforeFloor { node_id: next });
            }

            path.push(next);
            previous = Some(current);
            current = next;
        }
    }

    fn first_non_backtracking_outgoing_neighbor(
        &self,
        current: IsogenyGraphNodeId,
        previous: Option<IsogenyGraphNodeId>,
    ) -> Option<IsogenyGraphNodeId> {
        self.outgoing_edges(current)
            .map(|edge| edge.target())
            .find(|target| Some(*target) != previous)
    }
}
