use std::fmt;

use crate::fields::FieldError;

/// Errors produced by the educational `polynomials` module.
///
/// This enum centralizes the public failure modes that can arise in the
/// current polynomial APIs. Keeping them in one place avoids scattering raw
/// string literals through arithmetic, evaluation, interpolation, and
/// visualization helpers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PolynomialError {
    /// The chosen base field is structurally invalid for the requested
    /// polynomial algorithm.
    InvalidBaseField(FieldError),
    /// The current crate does not yet implement irreducibility testing for
    /// the requested base-field family.
    UnsupportedIrreducibilityBackend(&'static str),
    /// The current backend supports only a partial irreducibility procedure
    /// and the available exact criteria did not settle the input.
    UndeterminedIrreducibility(&'static str),
    /// Euclidean division was requested with the zero polynomial as divisor.
    DivisionByZeroPolynomial,
    /// A monic normalization was requested for the zero polynomial.
    ZeroPolynomialHasNoMonicNormalization,
    /// An algorithm expected a non-zero leading coefficient with an inverse
    /// but could not obtain one.
    NonInvertibleLeadingCoefficient,
    /// A multivariate term used an exponent vector whose arity does not match
    /// the declared polynomial arity.
    MonomialArityMismatch {
        /// The ambient arity declared for the polynomial.
        expected: usize,
        /// The actual arity carried by the monomial.
        actual: usize,
    },
    /// Two multivariate polynomials were combined even though they live in
    /// different ambient arities.
    IncompatibleMultivariateArity {
        /// Left-hand polynomial arity.
        lhs: usize,
        /// Right-hand polynomial arity.
        rhs: usize,
        /// Operation being attempted, such as `"addition"` or
        /// `"multiplication"`.
        operation: &'static str,
    },
    /// A multivariate evaluation point does not match the ambient arity of the
    /// polynomial being evaluated.
    EvaluationPointArityMismatch {
        /// Expected number of coordinates.
        expected: usize,
        /// Actual number of coordinates received.
        actual: usize,
    },
    /// Lagrange interpolation received repeated `x`-coordinates.
    DuplicateInterpolationAbscissa,
    /// Lagrange interpolation encountered a denominator that could not be
    /// inverted in the base field.
    NonInvertibleInterpolationDenominator,
}

impl fmt::Display for PolynomialError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBaseField(error) => {
                write!(
                    formatter,
                    "invalid base field for polynomial algorithm: {error}"
                )
            }
            Self::UnsupportedIrreducibilityBackend(message) => {
                write!(
                    formatter,
                    "irreducibility testing is not implemented for this base field: {message}"
                )
            }
            Self::UndeterminedIrreducibility(message) => write!(
                formatter,
                "irreducibility status could not be determined by the current exact partial backend: {message}"
            ),
            Self::DivisionByZeroPolynomial => {
                write!(formatter, "cannot divide by the zero polynomial")
            }
            Self::ZeroPolynomialHasNoMonicNormalization => {
                write!(formatter, "the zero polynomial has no monic normalization")
            }
            Self::NonInvertibleLeadingCoefficient => write!(
                formatter,
                "polynomial algorithm encountered a non-invertible leading coefficient"
            ),
            Self::MonomialArityMismatch { expected, actual } => write!(
                formatter,
                "all monomials must match the declared polynomial arity (expected {expected}, got {actual})"
            ),
            Self::IncompatibleMultivariateArity {
                lhs,
                rhs,
                operation,
            } => write!(
                formatter,
                "cannot perform multivariate {operation} with different arities (lhs {lhs}, rhs {rhs})"
            ),
            Self::EvaluationPointArityMismatch { expected, actual } => write!(
                formatter,
                "evaluation point arity must match the polynomial arity (expected {expected}, got {actual})"
            ),
            Self::DuplicateInterpolationAbscissa => write!(
                formatter,
                "lagrange interpolation requires distinct x-coordinates"
            ),
            Self::NonInvertibleInterpolationDenominator => write!(
                formatter,
                "lagrange interpolation encountered a non-invertible denominator"
            ),
        }
    }
}

impl std::error::Error for PolynomialError {}

#[cfg(test)]
mod tests {
    use super::PolynomialError;
    use crate::fields::FieldError;

    #[test]
    fn display_messages_remain_specific_to_the_polynomial_failure_mode() {
        assert_eq!(
            PolynomialError::InvalidBaseField(FieldError::InvalidModulus { modulus: 1 })
                .to_string(),
            "invalid base field for polynomial algorithm: invalid modulus: 1"
        );
        assert_eq!(
            PolynomialError::UnsupportedIrreducibilityBackend("complex approximate backend")
                .to_string(),
            "irreducibility testing is not implemented for this base field: complex approximate backend"
        );
        assert_eq!(
            PolynomialError::UndeterminedIrreducibility("no exact certificate found").to_string(),
            "irreducibility status could not be determined by the current exact partial backend: no exact certificate found"
        );
        assert_eq!(
            PolynomialError::DivisionByZeroPolynomial.to_string(),
            "cannot divide by the zero polynomial"
        );
        assert_eq!(
            PolynomialError::DuplicateInterpolationAbscissa.to_string(),
            "lagrange interpolation requires distinct x-coordinates"
        );
    }
}
