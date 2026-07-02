use core::fmt;

use num_bigint::BigUint;

use crate::elliptic_curves::CurveError;
use crate::fields::FieldError;

/// Placeholder error surface for future curve-isomorphism support.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CurveIsomorphismError {
    PointNotOnDomain,
    ImagePointNotOnCodomain,
    CurvesNotIsomorphic,
    NonInvertibleScale,
    UnsupportedCharacteristic { characteristic: BigUint },
    MissingSquareRoot,
    MissingFourthRoot,
    MissingSixthRoot,
    Field(FieldError),
    Curve(CurveError),
}

impl fmt::Display for CurveIsomorphismError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PointNotOnDomain => {
                write!(
                    formatter,
                    "the supplied point does not lie on the domain curve"
                )
            }
            Self::ImagePointNotOnCodomain => write!(
                formatter,
                "the computed image point does not lie on the codomain curve"
            ),
            Self::CurvesNotIsomorphic => write!(
                formatter,
                "the supplied short-Weierstrass curves are not isomorphic in the requested sense"
            ),
            Self::NonInvertibleScale => write!(
                formatter,
                "the proposed short-Weierstrass scaling factor is not invertible"
            ),
            Self::UnsupportedCharacteristic { characteristic } => write!(
                formatter,
                "short-Weierstrass isomorphism support requires characteristic different from 2 and 3, got {characteristic}"
            ),
            Self::MissingSquareRoot => write!(
                formatter,
                "the required square root does not exist in the current field"
            ),
            Self::MissingFourthRoot => write!(
                formatter,
                "the required fourth root does not exist in the current field"
            ),
            Self::MissingSixthRoot => write!(
                formatter,
                "the required sixth root does not exist in the current field"
            ),
            Self::Field(error) => write!(
                formatter,
                "field-level validation or arithmetic failed while working with a curve isomorphism: {error}"
            ),
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed while working with a curve isomorphism: {error}"
            ),
        }
    }
}

impl From<CurveError> for CurveIsomorphismError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl From<FieldError> for CurveIsomorphismError {
    fn from(error: FieldError) -> Self {
        Self::Field(error)
    }
}

impl std::error::Error for CurveIsomorphismError {}
