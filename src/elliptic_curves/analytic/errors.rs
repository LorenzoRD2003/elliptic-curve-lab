use core::fmt;

use crate::elliptic_curves::division_polynomials::DivisionPolynomialError;

/// Typed error surface for the educational complex-analytic elliptic-curve
/// milestone.
#[derive(Clone, Debug, PartialEq)]
pub enum AnalyticCurveError {
    TauNotInUpperHalfPlane,
    DegenerateLattice,
    NonPositiveLatticeOrientation,
    InvalidTorusTorsionIndex,
    InvalidEisensteinWeight,
    InvalidTruncationComparison,
    InvalidTruncationRadius,
    InvalidSeriesPrecision,
    NearlySingularAnalyticCurve,
    PointTooCloseToLatticePoint,
    PointNotInFundamentalParallelogram,
    InvalidModularMatrix,
    NonPositiveImaginaryPartAfterModularAction,
    UnsupportedNormalization,
    NumericalComparisonFailed,
}

impl fmt::Display for AnalyticCurveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            Self::TauNotInUpperHalfPlane => "tau is not in the upper half-plane",
            Self::DegenerateLattice => "lattice basis is degenerate",
            Self::NonPositiveLatticeOrientation => {
                "lattice basis does not have positive orientation"
            }
            Self::InvalidTorusTorsionIndex => {
                "torus torsion index must satisfy n > 0 and 0 ≤ a, b < n"
            }
            Self::InvalidEisensteinWeight => {
                "Eisenstein weight is invalid for the requested construction"
            }
            Self::InvalidTruncationComparison => {
                "truncation comparison requires a strictly larger second radius"
            }
            Self::InvalidTruncationRadius => "truncation radius must be positive and finite",
            Self::InvalidSeriesPrecision => "series precision parameter is invalid",
            Self::NearlySingularAnalyticCurve => "analytic Weierstrass model is nearly singular",
            Self::PointTooCloseToLatticePoint => "point is too close to a lattice point or pole",
            Self::PointNotInFundamentalParallelogram => {
                "point is not in the chosen fundamental parallelogram"
            }
            Self::InvalidModularMatrix => "matrix is not a valid modular transformation",
            Self::NonPositiveImaginaryPartAfterModularAction => {
                "modular action produced a non-positive imaginary part"
            }
            Self::UnsupportedNormalization => "requested normalization is not supported",
            Self::NumericalComparisonFailed => "numerical comparison or reduction failed",
        };

        write!(f, "{message}")
    }
}

impl std::error::Error for AnalyticCurveError {}

impl From<DivisionPolynomialError> for AnalyticCurveError {
    fn from(error: DivisionPolynomialError) -> Self {
        match error {
            DivisionPolynomialError::ZeroIndex => Self::InvalidTorusTorsionIndex,
            DivisionPolynomialError::Curve(crate::elliptic_curves::CurveError::SingularCurve) => {
                Self::NearlySingularAnalyticCurve
            }
            _ => Self::NumericalComparisonFailed,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AnalyticCurveError;
    use crate::elliptic_curves::{CurveError, division_polynomials::DivisionPolynomialError};

    #[test]
    fn display_for_upper_half_plane_error_is_human_readable() {
        assert_eq!(
            AnalyticCurveError::TauNotInUpperHalfPlane.to_string(),
            "tau is not in the upper half-plane"
        );
    }

    #[test]
    fn display_distinguishes_degenerate_and_orientation_failures() {
        assert_eq!(
            AnalyticCurveError::DegenerateLattice.to_string(),
            "lattice basis is degenerate"
        );
        assert_eq!(
            AnalyticCurveError::NonPositiveLatticeOrientation.to_string(),
            "lattice basis does not have positive orientation"
        );
    }

    #[test]
    fn display_for_invalid_torus_torsion_index_is_human_readable() {
        assert_eq!(
            AnalyticCurveError::InvalidTorusTorsionIndex.to_string(),
            "torus torsion index must satisfy n > 0 and 0 ≤ a, b < n"
        );
    }

    #[test]
    fn display_mentions_poles_for_near_lattice_point_failures() {
        assert_eq!(
            AnalyticCurveError::PointTooCloseToLatticePoint.to_string(),
            "point is too close to a lattice point or pole"
        );
    }

    #[test]
    fn division_polynomial_zero_index_maps_to_invalid_torus_torsion_index() {
        assert_eq!(
            AnalyticCurveError::from(DivisionPolynomialError::ZeroIndex),
            AnalyticCurveError::InvalidTorusTorsionIndex
        );
    }

    #[test]
    fn singular_curve_division_polynomial_errors_map_to_nearly_singular_curve() {
        assert_eq!(
            AnalyticCurveError::from(DivisionPolynomialError::Curve(CurveError::SingularCurve)),
            AnalyticCurveError::NearlySingularAnalyticCurve
        );
    }
}
