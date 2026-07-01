use core::fmt;

use crate::numerics::cornacchia::CornacchiaError;

/// Failure modes for exact binary-quadratic-form representation helpers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum QuadraticFormError {
    /// The diagonal coefficient `d` in `x² + d y²` must be positive.
    ZeroDiagonalCoefficient,
    /// The represented target must be positive for the current exact route.
    ZeroTarget,
    /// The current primitive-representation route failed in Cornacchia's
    /// algorithm.
    Cornacchia(CornacchiaError),
}

impl From<CornacchiaError> for QuadraticFormError {
    fn from(error: CornacchiaError) -> Self {
        match error {
            CornacchiaError::ZeroCoefficient => Self::ZeroDiagonalCoefficient,
            CornacchiaError::ZeroTarget | CornacchiaError::TrivialTarget => Self::ZeroTarget,
            other => Self::Cornacchia(other),
        }
    }
}

impl fmt::Display for QuadraticFormError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDiagonalCoefficient => write!(
                formatter,
                "the diagonal form x² + d y² requires a positive coefficient d"
            ),
            Self::ZeroTarget => write!(
                formatter,
                "quadratic-form representation currently requires a positive target"
            ),
            Self::Cornacchia(error) => write!(
                formatter,
                "Cornacchia-backed primitive representation failed: {error}"
            ),
        }
    }
}

impl std::error::Error for QuadraticFormError {}
