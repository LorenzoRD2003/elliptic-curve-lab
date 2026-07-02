use core::fmt;

use num_bigint::BigUint;

use crate::elliptic_curves::{CurveError, traits::CurveModel};
use crate::fields::FieldError;

/// Errors returned by explicit conversions between concrete curve models.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CurveModelConversionError {
    PointNotOnSource,
    PointNotOnTarget,
    UnsupportedCharacteristic { characteristic: BigUint },
    Field(FieldError),
    Curve(CurveError),
}

impl fmt::Display for CurveModelConversionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PointNotOnSource => {
                write!(
                    formatter,
                    "the supplied point does not lie on the source curve model"
                )
            }
            Self::PointNotOnTarget => {
                write!(
                    formatter,
                    "the supplied point does not lie on the target curve model"
                )
            }
            Self::UnsupportedCharacteristic { characteristic } => write!(
                formatter,
                "the requested curve-model conversion is not supported in characteristic {characteristic}"
            ),
            Self::Field(error) => write!(
                formatter,
                "field-level validation or arithmetic failed during the curve-model conversion: {error}"
            ),
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed during the curve-model conversion: {error}"
            ),
        }
    }
}

impl From<CurveError> for CurveModelConversionError {
    fn from(error: CurveError) -> Self {
        match error {
            CurveError::UnsupportedCharacteristic { characteristic } => {
                Self::UnsupportedCharacteristic { characteristic }
            }
            CurveError::Field(error) => Self::Field(error),
            other => Self::Curve(other),
        }
    }
}

impl From<FieldError> for CurveModelConversionError {
    fn from(error: FieldError) -> Self {
        Self::Field(error)
    }
}

impl From<CurveModelConversionError> for CurveError {
    fn from(error: CurveModelConversionError) -> Self {
        match error {
            CurveModelConversionError::PointNotOnSource
            | CurveModelConversionError::PointNotOnTarget => CurveError::PointNotOnCurve,
            CurveModelConversionError::UnsupportedCharacteristic { characteristic } => {
                CurveError::UnsupportedCharacteristic { characteristic }
            }
            CurveModelConversionError::Field(error) => CurveError::Field(error),
            CurveModelConversionError::Curve(error) => error,
        }
    }
}

impl std::error::Error for CurveModelConversionError {}

/// Explicit reusable witness for converting one concrete curve model into
/// another together with point transport in both directions.
pub trait CurveModelConversion {
    type Source: CurveModel;
    type Target: CurveModel;

    /// Returns the source curve model.
    fn source(&self) -> &Self::Source;

    /// Returns the target curve model.
    fn target(&self) -> &Self::Target;

    /// Transports one point from the source model to the target model.
    fn map_source_point(
        &self,
        point: &<Self::Source as CurveModel>::Point,
    ) -> Result<<Self::Target as CurveModel>::Point, CurveModelConversionError>;

    /// Transports one point from the target model back to the source model.
    fn map_target_point(
        &self,
        point: &<Self::Target as CurveModel>::Point,
    ) -> Result<<Self::Source as CurveModel>::Point, CurveModelConversionError>;
}

pub(crate) struct ReversedCurveModelConversion<C>(pub(crate) C);

impl<C> CurveModelConversion for ReversedCurveModelConversion<C>
where
    C: CurveModelConversion,
{
    type Source = C::Target;
    type Target = C::Source;

    fn source(&self) -> &Self::Source {
        self.0.target()
    }

    fn target(&self) -> &Self::Target {
        self.0.source()
    }

    fn map_source_point(
        &self,
        point: &<Self::Source as CurveModel>::Point,
    ) -> Result<<Self::Target as CurveModel>::Point, CurveModelConversionError> {
        self.0.map_target_point(point)
    }

    fn map_target_point(
        &self,
        point: &<Self::Target as CurveModel>::Point,
    ) -> Result<<Self::Source as CurveModel>::Point, CurveModelConversionError> {
        self.0.map_source_point(point)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{CurveError, traits::CurveModelConversionError};
    use num_bigint::BigUint;

    #[test]
    fn curve_error_conversion_preserves_unsupported_characteristic_as_a_first_class_variant() {
        let error = CurveModelConversionError::from(CurveError::UnsupportedCharacteristic {
            characteristic: BigUint::from(3u8),
        });

        assert_eq!(
            error,
            CurveModelConversionError::UnsupportedCharacteristic {
                characteristic: BigUint::from(3u8)
            }
        );
    }

    #[test]
    fn conversion_error_can_degrade_source_and_target_point_failures_to_point_not_on_curve() {
        assert_eq!(
            CurveError::from(CurveModelConversionError::PointNotOnSource),
            CurveError::PointNotOnCurve
        );
        assert_eq!(
            CurveError::from(CurveModelConversionError::PointNotOnTarget),
            CurveError::PointNotOnCurve
        );
    }
}
