use std::collections::HashSet;
use std::fmt;
use std::hash::Hash;

use num_bigint::BigUint;
use num_traits::One;

use crate::elliptic_curves::traits::PointIndexSampler;
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
    /// first stored outgoing edge that does not immediately backtrack.
    /// Use [`Self::find_floor_path_with_sampler`] for the randomized neighbor
    /// choice from Sutherland's `FindFloor`.
    ///
    /// This is not `FindShortestPathToFloor`, so
    /// [`FloorPathReport::distance_to_floor`] records only the path length
    /// found by this deterministic walk.
    pub fn find_floor_path(
        &self,
        start: IsogenyGraphNodeId,
        prime: &BigUint,
    ) -> Result<FloorPathReport, VolcanoSearchError> {
        self.find_floor_path_by_next_step(start, prime, |graph, current, previous| {
            graph
                .first_non_backtracking_outgoing_neighbor(current, previous)
                .ok_or(VolcanoSearchError::NoNonBacktrackingNeighbor { node_id: current })
        })
    }

    /// Follows Sutherland's randomized `FindFloor` neighbor-choice rule.
    ///
    /// The caller supplies a small index sampler instead of the crate depending
    /// on a randomness library. At the start vertex this samples any outgoing
    /// neighbor. At later vertices it samples uniformly from the stored
    /// outgoing neighbors except the immediate predecessor, matching the
    /// non-backtracking path convention in §3.1.
    ///
    /// As with [`Self::find_floor_path`], the result is a path found by this
    /// walk, not a certified shortest path to the floor.
    pub fn find_floor_path_with_sampler<S: PointIndexSampler>(
        &self,
        start: IsogenyGraphNodeId,
        prime: &BigUint,
        sampler: &mut S,
    ) -> Result<FloorPathReport, VolcanoSearchError> {
        self.find_floor_path_by_next_step(start, prime, |graph, current, previous| {
            graph.sample_non_backtracking_outgoing_neighbor(current, previous, sampler)
        })
    }

    fn find_floor_path_by_next_step<F>(
        &self,
        start: IsogenyGraphNodeId,
        prime: &BigUint,
        mut choose_next: F,
    ) -> Result<FloorPathReport, VolcanoSearchError>
    where
        F: FnMut(
            &Self,
            IsogenyGraphNodeId,
            Option<IsogenyGraphNodeId>,
        ) -> Result<IsogenyGraphNodeId, VolcanoSearchError>,
    {
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

            let next = choose_next(self, current, previous)?;

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

    fn sample_non_backtracking_outgoing_neighbor<S: PointIndexSampler>(
        &self,
        current: IsogenyGraphNodeId,
        previous: Option<IsogenyGraphNodeId>,
        sampler: &mut S,
    ) -> Result<IsogenyGraphNodeId, VolcanoSearchError> {
        let candidates = self
            .outgoing_edges(current)
            .map(|edge| edge.target())
            .filter(|target| Some(*target) != previous)
            .collect::<Vec<_>>();

        if candidates.is_empty() {
            return Err(VolcanoSearchError::NoNonBacktrackingNeighbor { node_id: current });
        }

        let upper_bound = candidates.len();
        let sampled_index = sampler
            .sample_index(upper_bound)
            .ok_or(VolcanoSearchError::SamplerExhausted { node_id: current })?;

        candidates
            .get(sampled_index)
            .copied()
            .ok_or(VolcanoSearchError::SamplerIndexOutOfRange {
                node_id: current,
                sampled_index,
                upper_bound,
            })
    }
}
