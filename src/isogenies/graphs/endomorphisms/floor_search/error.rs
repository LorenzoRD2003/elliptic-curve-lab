use std::fmt;

use num_bigint::BigUint;

use crate::isogenies::graphs::IsogenyGraphNodeId;
use crate::numerics::PositivePrimeError;

/// Errors produced while searching for floor evidence in an `ℓ`-volcano.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VolcanoSearchError {
    /// The chosen local parameter is not a positive prime `ℓ`.
    InvalidLocalPrime(PositivePrimeError),
    /// The requested node is not present in the stored graph.
    NodeNotFound { node_id: IsogenyGraphNodeId },
    /// The node was discovered, but its outgoing `ℓ`-isogenies were not fully
    /// expanded by the graph builder.
    NodeNotFullyExpanded { node_id: IsogenyGraphNodeId },
    /// The clean ordinary-volcano criterion excludes `j = 0` and `j = 1728`.
    SpecialJInvariant { node_id: IsogenyGraphNodeId },
    /// The complete local degree does not match either floor behavior
    /// `deg(v) ≤ 2` or non-floor behavior `deg(v) = ℓ + 1`.
    InconsistentWithVolcanoModel {
        node_id: IsogenyGraphNodeId,
        observed_out_degree: usize,
        expected_non_floor_degree: BigUint,
    },
    /// The current non-backtracking walk had no legal next edge.
    NoNonBacktrackingNeighbor { node_id: IsogenyGraphNodeId },
    /// The sampler did not provide an index for the available outgoing choices.
    SamplerExhausted { node_id: IsogenyGraphNodeId },
    /// The sampler returned an index outside `0..upper_bound`.
    SamplerIndexOutOfRange {
        node_id: IsogenyGraphNodeId,
        sampled_index: usize,
        upper_bound: usize,
    },
    /// The non-backtracking walk entered a cycle before seeing floor evidence.
    CycleDetectedBeforeFloor { node_id: IsogenyGraphNodeId },
    /// The shortest-path search exhausted its tracked candidate paths without
    /// certifying a floor vertex.
    NoFloorPathFound { start: IsogenyGraphNodeId },
}

impl fmt::Display for VolcanoSearchError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLocalPrime(error) => {
                write!(formatter, "invalid volcano search prime: {error}")
            }
            Self::NodeNotFound { node_id } => {
                write!(formatter, "node {:?} is not present in the graph", node_id)
            }
            Self::NodeNotFullyExpanded { node_id } => write!(
                formatter,
                "node {:?} was not fully expanded, so its observed outgoing degree is partial",
                node_id
            ),
            Self::SpecialJInvariant { node_id } => write!(
                formatter,
                "node {:?} has special j-invariant 0 or 1728, outside the clean ordinary-volcano criterion",
                node_id
            ),
            Self::InconsistentWithVolcanoModel {
                node_id,
                observed_out_degree,
                expected_non_floor_degree,
            } => write!(
                formatter,
                "node {:?} has complete outgoing degree {observed_out_degree}, neither floor-like nor equal to ℓ + 1 = {expected_non_floor_degree}",
                node_id
            ),
            Self::NoNonBacktrackingNeighbor { node_id } => write!(
                formatter,
                "node {:?} has no non-backtracking outgoing neighbor",
                node_id
            ),
            Self::SamplerExhausted { node_id } => write!(
                formatter,
                "the floor-search sampler did not choose an outgoing neighbor for node {:?}",
                node_id
            ),
            Self::SamplerIndexOutOfRange {
                node_id,
                sampled_index,
                upper_bound,
            } => write!(
                formatter,
                "the floor-search sampler returned index {sampled_index} for node {:?}, but the valid range is 0..{upper_bound}",
                node_id
            ),
            Self::CycleDetectedBeforeFloor { node_id } => write!(
                formatter,
                "floor search revisited node {:?} before finding floor evidence",
                node_id
            ),
            Self::NoFloorPathFound { start } => write!(
                formatter,
                "shortest floor search exhausted its tracked paths from {:?} without certifying a floor vertex",
                start
            ),
        }
    }
}

impl std::error::Error for VolcanoSearchError {}

impl From<PositivePrimeError> for VolcanoSearchError {
    fn from(error: PositivePrimeError) -> Self {
        Self::InvalidLocalPrime(error)
    }
}
