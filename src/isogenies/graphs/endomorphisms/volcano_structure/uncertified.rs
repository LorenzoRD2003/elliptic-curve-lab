use num_bigint::BigUint;

use crate::isogenies::graphs::{IsogenyGraphNodeId, endomorphisms::VolcanoSearchError};

/// Stored node whose volcano level could not be certified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UncertifiedVolcanoNodeReport {
    node_id: IsogenyGraphNodeId,
    reason: VolcanoStructureUncertifiedReason,
}

impl UncertifiedVolcanoNodeReport {
    pub(crate) fn new(
        node_id: IsogenyGraphNodeId,
        reason: VolcanoStructureUncertifiedReason,
    ) -> Self {
        Self { node_id, reason }
    }

    /// Returns the node that could not be placed on a certified level.
    pub fn node_id(&self) -> IsogenyGraphNodeId {
        self.node_id
    }

    /// Returns why the node was left outside the certified level partition.
    pub fn reason(&self) -> &VolcanoStructureUncertifiedReason {
        &self.reason
    }
}

/// Reason a stored graph node was not assigned a certified volcano level.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VolcanoStructureUncertifiedReason {
    /// The node's outgoing `ℓ`-isogenies were not fully expanded.
    PartialGraph,
    /// The node has `j = 0` or `j = 1728`, outside the clean criterion.
    SpecialJInvariant,
    /// The complete local degree does not match the ordinary-volcano model.
    InconsistentWithVolcanoModel {
        /// Complete outgoing degree observed in the stored graph.
        observed_out_degree: usize,
        /// Expected non-floor degree `ℓ + 1`.
        expected_non_floor_degree: BigUint,
    },
    /// The graph gave no legal non-backtracking next step toward the floor.
    NoNonBacktrackingNeighbor,
    /// A deterministic path search cycled before seeing floor evidence.
    CycleDetectedBeforeFloor,
    /// The shortest-path search exhausted its tracked paths.
    NoFloorPathFound,
}

impl VolcanoStructureUncertifiedReason {
    pub(crate) fn from_search_error(error: VolcanoSearchError) -> Result<Self, VolcanoSearchError> {
        match error {
            VolcanoSearchError::InvalidLocalPrime(error) => {
                Err(VolcanoSearchError::InvalidLocalPrime(error))
            }
            VolcanoSearchError::NodeNotFound { node_id } => {
                Err(VolcanoSearchError::NodeNotFound { node_id })
            }
            VolcanoSearchError::NodeNotFullyExpanded { .. } => Ok(Self::PartialGraph),
            VolcanoSearchError::SpecialJInvariant { .. } => Ok(Self::SpecialJInvariant),
            VolcanoSearchError::InconsistentWithVolcanoModel {
                observed_out_degree,
                expected_non_floor_degree,
                ..
            } => Ok(Self::InconsistentWithVolcanoModel {
                observed_out_degree,
                expected_non_floor_degree,
            }),
            VolcanoSearchError::NoNonBacktrackingNeighbor { .. } => {
                Ok(Self::NoNonBacktrackingNeighbor)
            }
            VolcanoSearchError::SamplerExhausted { node_id } => {
                Err(VolcanoSearchError::SamplerExhausted { node_id })
            }
            VolcanoSearchError::SamplerIndexOutOfRange {
                node_id,
                sampled_index,
                upper_bound,
            } => Err(VolcanoSearchError::SamplerIndexOutOfRange {
                node_id,
                sampled_index,
                upper_bound,
            }),
            VolcanoSearchError::CycleDetectedBeforeFloor { .. } => {
                Ok(Self::CycleDetectedBeforeFloor)
            }
            VolcanoSearchError::NoFloorPathFound { .. } => Ok(Self::NoFloorPathFound),
        }
    }
}
