use core::fmt;

use crate::isogenies::graphs::IsogenyGraphNodeId;

/// Failure modes for a user-supplied crater orientation witness.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CraterOrientationWitnessError {
    /// The crater report had no certified crater node to orient.
    EmptyCrater,
    /// A certified crater node had no declared successor.
    MissingSuccessor { source: IsogenyGraphNodeId },
    /// The declaration included a source outside the certified crater.
    SourceOutsideCrater { source: IsogenyGraphNodeId },
    /// A declared successor left the certified crater.
    TargetOutsideCrater {
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
    },
    /// The declared successor did not follow a certified internal horizontal edge.
    MissingCertifiedHorizontalEdge {
        source: IsogenyGraphNodeId,
        target: IsogenyGraphNodeId,
    },
    /// The declared successors did not close a single crater cycle from `start`.
    DoesNotCloseCycle { start: IsogenyGraphNodeId },
    /// The declared walk repeated a node before closing the crater cycle.
    RepeatsBeforeClosing { node: IsogenyGraphNodeId },
}

impl fmt::Display for CraterOrientationWitnessError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyCrater => write!(formatter, "cannot orient an empty certified crater"),
            Self::MissingSuccessor { source } => {
                write!(
                    formatter,
                    "node {:?} has no declared orientation successor",
                    source
                )
            }
            Self::SourceOutsideCrater { source } => {
                write!(
                    formatter,
                    "declared source {:?} is outside the certified crater",
                    source
                )
            }
            Self::TargetOutsideCrater { source, target } => write!(
                formatter,
                "declared successor {:?} for source {:?} is outside the certified crater",
                target, source
            ),
            Self::MissingCertifiedHorizontalEdge { source, target } => write!(
                formatter,
                "declared step {:?} -> {:?} is not a certified internal horizontal crater edge",
                source, target
            ),
            Self::DoesNotCloseCycle { start } => write!(
                formatter,
                "declared orientation does not close a crater cycle from {:?}",
                start
            ),
            Self::RepeatsBeforeClosing { node } => write!(
                formatter,
                "declared orientation repeats {:?} before closing the crater cycle",
                node
            ),
        }
    }
}

impl std::error::Error for CraterOrientationWitnessError {}
