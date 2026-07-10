use core::fmt;

/// Failure modes for binary-quadratic-form operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinaryQuadraticFormError {
    /// The requested operation requires a positive-definite form.
    NotPositiveDefinite,
    /// The form does not have the discriminant `D` of the quadratic class group.
    ///
    /// Class-group operations compose classes inside one fixed imaginary
    /// quadratic order, so every input representative must have that same
    /// discriminant.
    ClassGroupDiscriminantMismatch,
    /// The form is not primitive.
    ///
    /// The staged class-group layer works with primitive binary quadratic
    /// forms, equivalently proper invertible ideal classes for the order.
    NotPrimitive,
    /// The form is not the reduced positive-definite representative expected by the group API.
    ///
    /// Composition will normalize results by reduction, but the first public
    /// entry points use reduced representatives as their input convention.
    NotReducedPositiveDefinite,
    /// The first class-group enumerator requires a negative discriminant.
    NotNegativeDiscriminant,
    /// A quadratic-order discriminant must satisfy `D ≡ 0, 1 (mod 4)`.
    NotQuadraticOrderDiscriminant,
}

impl fmt::Display for BinaryQuadraticFormError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotPositiveDefinite => {
                write!(f, "binary quadratic form is not positive definite")
            }
            Self::ClassGroupDiscriminantMismatch => {
                write!(
                    f,
                    "binary quadratic form does not have the class-group discriminant"
                )
            }
            Self::NotPrimitive => write!(f, "binary quadratic form is not primitive"),
            Self::NotReducedPositiveDefinite => {
                write!(
                    f,
                    "binary quadratic form is not a reduced positive-definite representative"
                )
            }
            Self::NotNegativeDiscriminant => {
                write!(f, "quadratic class group enumeration requires D < 0")
            }
            Self::NotQuadraticOrderDiscriminant => {
                write!(
                    f,
                    "quadratic discriminant must be congruent to 0 or 1 modulo 4"
                )
            }
        }
    }
}

impl std::error::Error for BinaryQuadraticFormError {}
