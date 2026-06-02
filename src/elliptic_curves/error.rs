use core::fmt;

/// Errors returned when validating elliptic-curve models.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CurveError {
    /// The short Weierstrass model requires characteristic different from 2 and 3.
    UnsupportedCharacteristic { characteristic: u64 },
    /// A torsion helper received an invalid order parameter.
    InvalidTorsionOrder { order: usize },
    /// The supplied coefficients define a singular cubic.
    SingularCurve,
    /// The supplied affine coordinates do not satisfy the curve equation.
    PointNotOnCurve,
    /// An exhaustively checked finite-group axiom failed.
    GroupAxiomViolation { axiom: &'static str },
}

impl fmt::Display for CurveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedCharacteristic { characteristic } => write!(
                f,
                "short Weierstrass form requires characteristic different from 2 and 3, got {characteristic}"
            ),
            Self::InvalidTorsionOrder { order } => {
                write!(f, "torsion order must be a positive integer, got {order}")
            }
            Self::SingularCurve => {
                write!(f, "short Weierstrass coefficients define a singular curve")
            }
            Self::PointNotOnCurve => {
                write!(f, "affine coordinates do not satisfy the curve equation")
            }
            Self::GroupAxiomViolation { axiom } => {
                write!(f, "finite-group axiom validation failed: {axiom}")
            }
        }
    }
}

impl std::error::Error for CurveError {}
