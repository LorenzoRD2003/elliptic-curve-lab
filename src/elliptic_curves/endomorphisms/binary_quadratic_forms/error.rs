use core::fmt;

/// Failure modes for binary-quadratic-form operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryQuadraticFormError {
    /// The requested operation requires a positive-definite form.
    NotPositiveDefinite,
}

impl fmt::Display for BinaryQuadraticFormError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPositiveDefinite => {
                write!(f, "binary quadratic form is not positive definite")
            }
        }
    }
}

impl std::error::Error for BinaryQuadraticFormError {}
