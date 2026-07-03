use core::fmt;

use crate::elliptic_curves::{CurveError, short_weierstrass::isomorphisms::CurveIsomorphismError};

/// Errors produced by the staged rational-torsion workflow over `Q`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum RationalTorsionError {
    /// A curve-level operation failed while preparing or checking candidates.
    Curve(CurveError),
    /// A scaling-isomorphism operation failed while transporting to an integral model.
    Isomorphism(CurveIsomorphismError),
    /// The current implementation cannot yet produce a certified integral model.
    IntegralModelUnavailable,
    /// A requested group shape is not allowed by Mazur's theorem over `Q`.
    InvalidMazurShape {
        /// The requested Mazur family.
        family: &'static str,
        /// The invalid family parameter.
        value: usize,
    },
    /// Candidate verification produced a torsion shape outside Mazur's theorem.
    InconsistentMazurShape {
        /// Number of points found before the impossible shape was detected.
        point_count: usize,
    },
}

impl fmt::Display for RationalTorsionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed while working with rational torsion: {error}"
            ),
            Self::Isomorphism(error) => write!(
                formatter,
                "short-Weierstrass scaling failed while working with rational torsion: {error}"
            ),
            Self::IntegralModelUnavailable => write!(
                formatter,
                "the current rational-torsion route cannot yet certify an integral short-Weierstrass model"
            ),
            Self::InvalidMazurShape { family, value } => write!(
                formatter,
                "the requested {family} rational torsion shape with parameter {value} is not allowed by Mazur's theorem over Q"
            ),
            Self::InconsistentMazurShape { point_count } => write!(
                formatter,
                "candidate verification found {point_count} torsion points, which does not match a Mazur torsion shape"
            ),
        }
    }
}

impl From<CurveError> for RationalTorsionError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl From<CurveIsomorphismError> for RationalTorsionError {
    fn from(error: CurveIsomorphismError) -> Self {
        Self::Isomorphism(error)
    }
}

impl std::error::Error for RationalTorsionError {}
