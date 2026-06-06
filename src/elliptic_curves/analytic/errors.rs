use core::fmt;

use crate::elliptic_curves::division_polynomials::DivisionPolynomialError;
use crate::numerics::SimpsonIntegrationError;

/// Typed error surface for the educational complex-analytic elliptic-curve
/// layer.
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
    CubicRootRecoveryFailed,
    RepeatedCubicRoot,
    AmbiguousRootOrdering,
    InvalidLegendreModulus,
    InvalidAgmInput,
    InvalidEllipticIntegralInput,
    InvalidPeriodRecoveryConfig,
    PeriodRecoveryFailed,
    PeriodRatioNotInUpperHalfPlane,
    InverseUniformizationFailed,
    AbelJacobiIntegrationFailed,
    BranchChoiceAmbiguous,
    PeriodValidationFailed,
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
            Self::CubicRootRecoveryFailed => "failed to recover the cubic roots stably",
            Self::RepeatedCubicRoot => "cubic root recovery produced a repeated root",
            Self::AmbiguousRootOrdering => {
                "cubic roots could not be ordered consistently for the requested normalization"
            }
            Self::InvalidLegendreModulus => {
                "Legendre modulus is invalid for the requested period construction"
            }
            Self::InvalidAgmInput => "AGM input is invalid for the requested numerical iteration",
            Self::InvalidEllipticIntegralInput => {
                "elliptic-integral input is invalid for the requested branch or domain"
            }
            Self::InvalidPeriodRecoveryConfig => "period-recovery configuration is invalid",
            Self::PeriodRecoveryFailed => "period recovery failed",
            Self::PeriodRatioNotInUpperHalfPlane => {
                "recovered period ratio does not lie in the upper half-plane"
            }
            Self::InverseUniformizationFailed => "inverse uniformization failed",
            Self::AbelJacobiIntegrationFailed => "Abel-Jacobi integration failed",
            Self::BranchChoiceAmbiguous => {
                "branch choice is ambiguous for the requested analytic continuation"
            }
            Self::PeriodValidationFailed => {
                "recovered periods failed the requested validation checks"
            }
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

impl From<SimpsonIntegrationError<AnalyticCurveError>> for AnalyticCurveError {
    fn from(error: SimpsonIntegrationError<AnalyticCurveError>) -> Self {
        match error {
            SimpsonIntegrationError::NonFiniteIntegrandValue { .. } => {
                Self::AbelJacobiIntegrationFailed
            }
            SimpsonIntegrationError::Integrand(error) => error,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AnalyticCurveError;
    use crate::elliptic_curves::{CurveError, division_polynomials::DivisionPolynomialError};
    use crate::numerics::SimpsonIntegrationError;

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
    fn display_for_period_recovery_errors_is_mathematically_specific() {
        assert_eq!(
            AnalyticCurveError::CubicRootRecoveryFailed.to_string(),
            "failed to recover the cubic roots stably"
        );
        assert_eq!(
            AnalyticCurveError::RepeatedCubicRoot.to_string(),
            "cubic root recovery produced a repeated root"
        );
        assert_eq!(
            AnalyticCurveError::AmbiguousRootOrdering.to_string(),
            "cubic roots could not be ordered consistently for the requested normalization"
        );
        assert_eq!(
            AnalyticCurveError::InvalidLegendreModulus.to_string(),
            "Legendre modulus is invalid for the requested period construction"
        );
        assert_eq!(
            AnalyticCurveError::InvalidAgmInput.to_string(),
            "AGM input is invalid for the requested numerical iteration"
        );
        assert_eq!(
            AnalyticCurveError::InvalidEllipticIntegralInput.to_string(),
            "elliptic-integral input is invalid for the requested branch or domain"
        );
        assert_eq!(
            AnalyticCurveError::InvalidPeriodRecoveryConfig.to_string(),
            "period-recovery configuration is invalid"
        );
        assert_eq!(
            AnalyticCurveError::PeriodRecoveryFailed.to_string(),
            "period recovery failed"
        );
        assert_eq!(
            AnalyticCurveError::PeriodRatioNotInUpperHalfPlane.to_string(),
            "recovered period ratio does not lie in the upper half-plane"
        );
        assert_eq!(
            AnalyticCurveError::InverseUniformizationFailed.to_string(),
            "inverse uniformization failed"
        );
        assert_eq!(
            AnalyticCurveError::AbelJacobiIntegrationFailed.to_string(),
            "Abel-Jacobi integration failed"
        );
        assert_eq!(
            AnalyticCurveError::BranchChoiceAmbiguous.to_string(),
            "branch choice is ambiguous for the requested analytic continuation"
        );
        assert_eq!(
            AnalyticCurveError::PeriodValidationFailed.to_string(),
            "recovered periods failed the requested validation checks"
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

    #[test]
    fn simpson_non_finite_integrand_value_maps_to_abel_jacobi_integration_failure() {
        assert_eq!(
            AnalyticCurveError::from(
                SimpsonIntegrationError::<AnalyticCurveError>::NonFiniteIntegrandValue {
                    index: 3,
                    parameter: 0.25,
                }
            ),
            AnalyticCurveError::AbelJacobiIntegrationFailed
        );
    }

    #[test]
    fn simpson_integrand_error_maps_through_transparently() {
        assert_eq!(
            AnalyticCurveError::from(SimpsonIntegrationError::Integrand(
                AnalyticCurveError::BranchChoiceAmbiguous
            )),
            AnalyticCurveError::BranchChoiceAmbiguous
        );
    }
}
