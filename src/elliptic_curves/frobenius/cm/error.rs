use core::fmt;

use crate::{
    elliptic_curves::CurveError,
    numerics::{cornacchia::CornacchiaError, quadratic_forms::QuadraticFormError},
};

/// Failure modes for CM trace-candidate generation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CmTraceCandidateError {
    /// The CM discriminant must be negative.
    NonNegativeDiscriminant,
    /// The prime-like input `p` must be positive.
    ZeroPrime,
    /// The underlying Cornacchia candidate route failed.
    Cornacchia(CornacchiaError),
    /// The auxiliary primitive representation route for `p` failed.
    QuadraticForm(QuadraticFormError),
    /// A curve-side group-law operation failed while testing a trace sign.
    Curve(CurveError),
}

impl From<CornacchiaError> for CmTraceCandidateError {
    fn from(error: CornacchiaError) -> Self {
        Self::Cornacchia(error)
    }
}

impl From<QuadraticFormError> for CmTraceCandidateError {
    fn from(error: QuadraticFormError) -> Self {
        Self::QuadraticForm(error)
    }
}

impl From<CurveError> for CmTraceCandidateError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl fmt::Display for CmTraceCandidateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonNegativeDiscriminant => write!(
                formatter,
                "CM trace candidates require a negative discriminant D"
            ),
            Self::ZeroPrime => write!(
                formatter,
                "CM trace candidates require a positive prime-like input p"
            ),
            Self::Cornacchia(error) => {
                write!(formatter, "Cornacchia candidate generation failed: {error}")
            }
            Self::QuadraticForm(error) => {
                write!(formatter, "quadratic-form representation failed: {error}")
            }
            Self::Curve(error) => write!(formatter, "curve-side trace-sign test failed: {error}"),
        }
    }
}

impl std::error::Error for CmTraceCandidateError {}
