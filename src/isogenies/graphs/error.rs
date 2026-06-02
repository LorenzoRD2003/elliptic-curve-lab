use core::fmt;

use crate::elliptic_curves::CurveError;
use crate::isogenies::IsogenyError;
use crate::isogenies::graphs::IsogenyGraphNodeId;

/// Errors produced by the milestone-6 educational isogeny-graph layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyGraphError {
    /// The requested isogeny degree is not valid for the current operation.
    InvalidDegree,
    /// The current builder path only accepts prime degrees.
    DegreeMustBePrimeForThisBuilder { degree: usize },
    /// A curve-domain operation failed while building or checking the graph.
    Curve(CurveError),
    /// An isogeny-domain operation failed while building or checking the graph.
    Isogeny(IsogenyError),
    /// The requested source node id is missing from the graph.
    MissingSourceNode(IsogenyGraphNodeId),
    /// The requested target node id is missing from the graph.
    MissingTargetNode(IsogenyGraphNodeId),
    /// No rational kernel of the requested degree was found in the current
    /// small finite setting.
    NonRationalKernelForRequestedDegree { degree: usize },
}

impl fmt::Display for IsogenyGraphError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDegree => write!(formatter, "the requested isogeny degree is invalid"),
            Self::DegreeMustBePrimeForThisBuilder { degree } => write!(
                formatter,
                "the current educational graph builder only supports prime degree {degree}"
            ),
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed while building or querying the isogeny graph: {error}"
            ),
            Self::Isogeny(error) => write!(
                formatter,
                "isogeny construction or validation failed while building the graph: {error}"
            ),
            Self::MissingSourceNode(id) => {
                write!(formatter, "missing source graph node with id {}", id.0)
            }
            Self::MissingTargetNode(id) => {
                write!(formatter, "missing target graph node with id {}", id.0)
            }
            Self::NonRationalKernelForRequestedDegree { degree } => write!(
                formatter,
                "no rational kernel of degree {degree} was found in the current small finite setting"
            ),
        }
    }
}

impl From<CurveError> for IsogenyGraphError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl From<IsogenyError> for IsogenyGraphError {
    fn from(error: IsogenyError) -> Self {
        Self::Isogeny(error)
    }
}

impl std::error::Error for IsogenyGraphError {}

#[cfg(test)]
mod tests {
    use super::IsogenyGraphError;
    use crate::elliptic_curves::CurveError;
    use crate::isogenies::IsogenyError;
    use crate::isogenies::graphs::IsogenyGraphNodeId;

    #[test]
    fn converts_curve_errors_into_graph_errors() {
        let error = IsogenyGraphError::from(CurveError::PointNotOnCurve);

        assert_eq!(error, IsogenyGraphError::Curve(CurveError::PointNotOnCurve));
    }

    #[test]
    fn converts_isogeny_errors_into_graph_errors() {
        let error = IsogenyGraphError::from(IsogenyError::EmptyKernel);

        assert_eq!(error, IsogenyGraphError::Isogeny(IsogenyError::EmptyKernel));
    }

    #[test]
    fn missing_node_errors_keep_the_requested_id() {
        assert_eq!(
            IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(4)),
            IsogenyGraphError::MissingSourceNode(IsogenyGraphNodeId(4))
        );
        assert_eq!(
            IsogenyGraphError::MissingTargetNode(IsogenyGraphNodeId(7)),
            IsogenyGraphError::MissingTargetNode(IsogenyGraphNodeId(7))
        );
    }
}
