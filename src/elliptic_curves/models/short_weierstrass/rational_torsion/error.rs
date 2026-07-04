use core::fmt;

use crate::elliptic_curves::{CurveError, short_weierstrass::isomorphisms::CurveIsomorphismError};

/// Errors produced by the staged rational-torsion workflow over `Q`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RationalTorsionError {
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
    /// A completed report paired a group with the wrong number of points.
    InconsistentReportGroup {
        /// Cardinality predicted by the classified torsion group.
        group_cardinality: usize,
        /// Number of points stored in the report.
        point_count: usize,
    },
    /// A completed report recorded fewer candidates than accepted points.
    InvalidCandidateAccounting {
        /// Number of candidates checked.
        candidate_count: usize,
        /// Number of accepted points stored in the report.
        point_count: usize,
    },
    /// The good-reduction/Hensel strategy could not complete.
    ReductionHenselUnavailable {
        /// Stable explanation of the failed internal stage.
        reason: &'static str,
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
            Self::InconsistentReportGroup {
                group_cardinality,
                point_count,
            } => write!(
                formatter,
                "rational-torsion report stores {point_count} points but the classified group has cardinality {group_cardinality}"
            ),
            Self::InvalidCandidateAccounting {
                candidate_count,
                point_count,
            } => write!(
                formatter,
                "rational-torsion report checked {candidate_count} candidates but stores {point_count} accepted points"
            ),
            Self::ReductionHenselUnavailable { reason } => write!(
                formatter,
                "the good-reduction/Hensel rational-torsion strategy could not complete: {reason}"
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
