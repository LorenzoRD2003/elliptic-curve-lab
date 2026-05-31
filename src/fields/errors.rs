use core::fmt;

/// Errors shared across field construction and arithmetic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FieldError {
    /// Division or inversion was requested with the additive identity.
    DivisionByZero,
    /// The configured modulus is invalid for the intended field family.
    InvalidModulus { modulus: u64 },
    /// The supplied polynomial modulus is structurally invalid.
    InvalidPolynomialModulus,
    /// The configured polynomial modulus is not irreducible.
    NonIrreduciblePolynomial,
    /// A value is outside the accepted canonical range.
    ElementOutOfRange { value: String },
    /// The operation mixes incompatible field descriptors.
    IncompatibleFieldParameters,
    /// The requested quantity does not fit the target representation.
    CardinalityOverflow,
    /// Placeholder for scaffolding paths that still need a specific error.
    Unsupported(&'static str),
}

impl fmt::Display for FieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DivisionByZero => write!(f, "division by zero"),
            Self::InvalidModulus { modulus } => write!(f, "invalid modulus: {modulus}"),
            Self::InvalidPolynomialModulus => write!(f, "invalid polynomial modulus"),
            Self::NonIrreduciblePolynomial => write!(f, "polynomial modulus is not irreducible"),
            Self::ElementOutOfRange { value } => write!(f, "element out of range: {value}"),
            Self::IncompatibleFieldParameters => {
                write!(f, "incompatible field parameters")
            }
            Self::CardinalityOverflow => write!(f, "field cardinality does not fit in u128"),
            Self::Unsupported(message) => write!(f, "unsupported operation: {message}"),
        }
    }
}

impl std::error::Error for FieldError {}
