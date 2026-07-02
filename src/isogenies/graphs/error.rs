use core::fmt;

use crate::elliptic_curves::{
    CurveError,
    endomorphisms::quadratic_orders::{
        ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorizationError,
    },
};
use crate::isogenies::{error::IsogenyError, graphs::IsogenyGraphNodeId};

/// Errors produced by the educational isogeny-graph layer.
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
    /// An endomorphism-side candidate-order computation failed for a stored node.
    Endomorphism(ImaginaryQuadraticOrderError),
    /// A Frobenius discriminant could not be factored into candidate quadratic-order data.
    EndomorphismFactorization(QuadraticDiscriminantFactorizationError),
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
            Self::Endomorphism(error) => write!(
                formatter,
                "endomorphism-side candidate-order derivation failed while querying the isogeny graph: {error}"
            ),
            Self::EndomorphismFactorization(error) => write!(
                formatter,
                "Frobenius discriminant factorization failed while deriving endomorphism-side candidates for the isogeny graph: {error}"
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

    use crate::elliptic_curves::{
        CurveError,
        endomorphisms::quadratic_orders::{
            ImaginaryQuadraticOrderError, QuadraticDiscriminantFactorizationError,
        },
    };
    use crate::isogenies::{
        error::{IsogenyError, IsogenyKernelError},
        graphs::{IsogenyGraphError, IsogenyGraphNodeId},
    };

    #[test]
    fn converts_curve_errors_into_graph_errors() {
        let error = IsogenyGraphError::from(CurveError::PointNotOnCurve);

        assert_eq!(error, IsogenyGraphError::Curve(CurveError::PointNotOnCurve));
    }

    #[test]
    fn converts_isogeny_errors_into_graph_errors() {
        let error = IsogenyGraphError::from(IsogenyError::Kernel(IsogenyKernelError::EmptyKernel));

        assert_eq!(
            error,
            IsogenyGraphError::Isogeny(IsogenyError::Kernel(IsogenyKernelError::EmptyKernel))
        );
    }

    #[test]
    fn converts_endomorphism_errors_into_graph_errors() {
        let error =
            IsogenyGraphError::from(ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant);

        assert_eq!(
            error,
            IsogenyGraphError::Endomorphism(
                ImaginaryQuadraticOrderError::NonImaginaryOrderDiscriminant
            )
        );
    }

    #[test]
    fn converts_endomorphism_factorization_errors_into_graph_errors() {
        let error =
            IsogenyGraphError::from(QuadraticDiscriminantFactorizationError::ZeroDiscriminant);

        assert_eq!(
            error,
            IsogenyGraphError::EndomorphismFactorization(
                QuadraticDiscriminantFactorizationError::ZeroDiscriminant
            )
        );
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
