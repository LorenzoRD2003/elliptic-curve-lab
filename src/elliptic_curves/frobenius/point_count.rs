use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::{CharacterSumPointCount, FrobeniusTrace, HasseInterval};

/// Public strategy choices for counting `#E(F_q)` on small finite curves.
///
/// The current educational implementation distinguishes:
///
/// - [`Self::Exhaustive`], which materializes all rational points directly
/// - [`Self::QuadraticCharacter`], which uses `#E(F_q) = q + 1 + Σ_{x ∈ F_q} χ(f(x))`
/// - [`Self::Auto`], which chooses the best implemented route for the curve
///
/// At the moment [`Self::Auto`] prefers the quadratic-character route when the
/// curve implements it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointCountStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
}

/// Shared point-count result returned by curve-side counting methods.
///
/// This enum keeps the algorithmic route explicit instead of collapsing every
/// count immediately into a bare integer. Callers can still recover the common
/// arithmetic data such as `#E(F_q)`, the trace `t`, and the Hasse interval
/// without losing which strategy produced it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PointCountReport {
    ExhaustiveTrace(FrobeniusTrace),
    QuadraticCharacter(CharacterSumPointCount),
}

impl PointCountReport {
    /// Returns the strategy used to build this report.
    pub fn strategy(&self) -> PointCountStrategy {
        match self {
            Self::ExhaustiveTrace(_) => PointCountStrategy::Exhaustive,
            Self::QuadraticCharacter(_) => PointCountStrategy::QuadraticCharacter,
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => trace.field_order(),
            Self::QuadraticCharacter(report) => report.field_order(),
        }
    }

    /// Returns the counted curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => u128::from(trace.curve_order()),
            Self::QuadraticCharacter(report) => report.curve_order(),
        }
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> i128 {
        match self {
            Self::ExhaustiveTrace(trace) => i128::from(trace.trace()),
            Self::QuadraticCharacter(report) => report.trace(),
        }
    }

    /// Returns the discrete Hasse interval attached to the same `F_q`.
    pub fn hasse_interval(&self) -> HasseInterval {
        match self {
            Self::ExhaustiveTrace(trace) => trace.hasse_interval(),
            Self::QuadraticCharacter(report) => report.hasse_interval(),
        }
    }

    /// Converts this report into the shared Frobenius-trace package.
    pub fn to_frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        match self {
            Self::ExhaustiveTrace(trace) => Ok(trace.clone()),
            Self::QuadraticCharacter(report) => report.to_frobenius_trace(),
        }
    }
}
