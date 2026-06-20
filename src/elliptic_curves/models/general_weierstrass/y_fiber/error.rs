use core::fmt;

use crate::elliptic_curves::CurveError;
use crate::fields::FieldError;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum GeneralWeierstrassYFiberError {
    UnsupportedCharacteristic { characteristic: u64 },
    Field(FieldError),
}

impl fmt::Display for GeneralWeierstrassYFiberError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedCharacteristic { characteristic } => write!(
                formatter,
                "solving y^2 + uy = v by completing the square is not supported in characteristic {characteristic}"
            ),
            Self::Field(error) => write!(
                formatter,
                "field-level arithmetic failed while preparing the y-quadratic solve: {error}"
            ),
        }
    }
}

impl std::error::Error for GeneralWeierstrassYFiberError {}

impl From<FieldError> for GeneralWeierstrassYFiberError {
    fn from(error: FieldError) -> Self {
        Self::Field(error)
    }
}

impl From<GeneralWeierstrassYFiberError> for CurveError {
    fn from(error: GeneralWeierstrassYFiberError) -> Self {
        match error {
            GeneralWeierstrassYFiberError::UnsupportedCharacteristic { characteristic } => {
                Self::UnsupportedCharacteristic { characteristic }
            }
            GeneralWeierstrassYFiberError::Field(error) => Self::Field(error),
        }
    }
}
