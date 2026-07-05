use std::fmt;

use crate::isogenies::graphs::{
    IsogenyGraphNodeId, endomorphisms::refinement::CandidateRefinementStrategy,
};

/// Errors produced while refining endomorphism candidates from graph evidence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CandidateRefinementError {
    /// The requested node is not present in the endomorphism report.
    NodeNotFound { node_id: IsogenyGraphNodeId },
    /// The chosen local parameter is not a positive prime `ℓ`.
    InvalidLocalPrime,
    /// The requested strategy is part of the staged API but has not been wired
    /// to executable refinement logic yet.
    StrategyNotImplemented {
        /// Strategy requested by the caller.
        strategy: CandidateRefinementStrategy,
    },
}

impl fmt::Display for CandidateRefinementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NodeNotFound { node_id } => {
                write!(
                    f,
                    "node {:?} is not present in the endomorphism report",
                    node_id
                )
            }
            Self::InvalidLocalPrime => write!(f, "the local parameter is not a positive prime ℓ"),
            Self::StrategyNotImplemented { strategy } => write!(
                f,
                "candidate refinement strategy {:?} is not implemented yet",
                strategy
            ),
        }
    }
}

impl std::error::Error for CandidateRefinementError {}
