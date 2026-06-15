use core::fmt;

use crate::elliptic_curves::CurveError;
use crate::polynomials::PolynomialError;

/// Errors produced by educational division-polynomial and torsion helpers for
/// short-Weierstrass curves.
///
/// This error surface is intentionally separate from `PolynomialError`
/// because several division-polynomial failure modes are not purely polynomial.
/// For example, asking for `ψ_2` as a polynomial only in `x` loses the
/// necessary `y`-factor and should therefore report a division-polynomial
/// modeling error rather than a generic polynomial arithmetic failure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DivisionPolynomialError {
    /// Division polynomials are indexed from `1` upward.
    ZeroIndex,
    /// The requested index is mathematically meaningful, but the current
    /// educational implementation does not yet support it.
    UnsupportedIndex {
        /// The requested division-polynomial index.
        n: usize,
    },
    /// Even division polynomials are not polynomials purely in the
    /// `x`-coordinate; they carry an extra `y`-factor.
    EvenIndexRequiresYFactor {
        /// The requested even division-polynomial index.
        n: usize,
    },
    /// A curve-level validation or geometric precondition failed.
    Curve(CurveError),
    /// The current educational evaluation helpers work only on finite affine
    /// points, not on the point at infinity.
    PointAtInfinityNotSupported,
    /// An underlying polynomial-domain operation failed.
    Polynomial(PolynomialError),
    /// The requested helper needs exhaustive access to all field elements,
    /// but the current base field does not expose honest enumeration.
    FieldNotEnumerable,
    /// The requested helper needs an honest square-root backend, but the
    /// current base field does not provide one.
    FieldHasNoSquareRootBackend,
}

impl fmt::Display for DivisionPolynomialError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroIndex => write!(
                formatter,
                "division polynomials are indexed from 1 upward; index 0 is undefined"
            ),
            Self::UnsupportedIndex { n } => write!(
                formatter,
                "division polynomial ψ_{n} is not supported by the current educational implementation"
            ),
            Self::EvenIndexRequiresYFactor { n } => write!(
                formatter,
                "division polynomial ψ_{n} is not a polynomial purely in x; the even-index formula requires a y-factor"
            ),
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed while working with a division polynomial or torsion helper: {error}"
            ),
            Self::PointAtInfinityNotSupported => write!(
                formatter,
                "division-polynomial evaluation at the point at infinity is not supported by the current affine educational API"
            ),
            Self::Polynomial(error) => write!(
                formatter,
                "polynomial operation failed while working with a division polynomial or torsion helper: {error}"
            ),
            Self::FieldNotEnumerable => write!(
                formatter,
                "this helper requires a base field that can enumerate all of its elements honestly"
            ),
            Self::FieldHasNoSquareRootBackend => write!(
                formatter,
                "this helper requires a base field with an honest square-root backend"
            ),
        }
    }
}

impl From<CurveError> for DivisionPolynomialError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl From<PolynomialError> for DivisionPolynomialError {
    fn from(error: PolynomialError) -> Self {
        Self::Polynomial(error)
    }
}

impl std::error::Error for DivisionPolynomialError {}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::CurveError;
    use crate::polynomials::PolynomialError;

    use super::DivisionPolynomialError;

    #[test]
    fn display_explains_even_index_x_only_limitation() {
        let error = DivisionPolynomialError::EvenIndexRequiresYFactor { n: 2 };
        let rendered = error.to_string();
        assert!(rendered.contains("ψ_2"));
        assert!(rendered.contains("requires a y-factor"));
    }

    #[test]
    fn from_curve_error_wraps_the_curve_variant() {
        let error = DivisionPolynomialError::from(CurveError::PointNotOnCurve);
        assert_eq!(
            error,
            DivisionPolynomialError::Curve(CurveError::PointNotOnCurve)
        );
    }

    #[test]
    fn from_polynomial_error_wraps_the_polynomial_variant() {
        let error = DivisionPolynomialError::from(PolynomialError::DivisionByZeroPolynomial);
        assert_eq!(
            error,
            DivisionPolynomialError::Polynomial(PolynomialError::DivisionByZeroPolynomial)
        );
    }
}
