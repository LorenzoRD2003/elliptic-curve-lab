use core::fmt;

/// Errors shared across field construction and arithmetic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FieldError {
    /// Division or inversion was requested with the additive identity.
    DivisionByZero,
    /// The configured modulus is invalid for the intended field family.
    InvalidModulus { modulus: String },
    /// The supplied polynomial modulus is structurally invalid.
    InvalidPolynomialModulus,
    /// The configured polynomial modulus is not irreducible.
    NonIrreduciblePolynomial,
    /// A non-zero quotient representative is not invertible.
    NonInvertibleElement,
    /// The current exact checks could not determine whether the modulus is
    /// irreducible.
    UndeterminedPolynomialModulusIrreducibility,
    /// A value is outside the accepted canonical range.
    ElementOutOfRange { value: String },
    /// The operation mixes incompatible field descriptors.
    IncompatibleFieldParameters,
    /// The requested quantity does not fit the target representation.
    CardinalityOverflow,
    /// The field characteristic does not fit the target representation.
    CharacteristicOverflow,
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
            Self::NonInvertibleElement => write!(f, "element is not invertible"),
            Self::UndeterminedPolynomialModulusIrreducibility => write!(
                f,
                "polynomial modulus irreducibility could not be determined by the current exact backend"
            ),
            Self::ElementOutOfRange { value } => write!(f, "element out of range: {value}"),
            Self::IncompatibleFieldParameters => {
                write!(f, "incompatible field parameters")
            }
            Self::CardinalityOverflow => {
                write!(
                    f,
                    "field cardinality does not fit the target representation"
                )
            }
            Self::CharacteristicOverflow => {
                write!(
                    f,
                    "field characteristic does not fit the target representation"
                )
            }
            Self::Unsupported(message) => write!(f, "unsupported operation: {message}"),
        }
    }
}

impl std::error::Error for FieldError {}

#[cfg(test)]
mod tests {

    use crate::fields::FieldError;

    #[test]
    fn display_messages_remain_specific_to_the_error_variant() {
        assert_eq!(FieldError::DivisionByZero.to_string(), "division by zero");
        assert_eq!(
            FieldError::InvalidModulus {
                modulus: "1".into()
            }
            .to_string(),
            "invalid modulus: 1"
        );
        assert_eq!(
            FieldError::ElementOutOfRange { value: "42".into() }.to_string(),
            "element out of range: 42"
        );
        assert_eq!(
            FieldError::Unsupported("cube roots over this backend").to_string(),
            "unsupported operation: cube roots over this backend"
        );
    }
}
