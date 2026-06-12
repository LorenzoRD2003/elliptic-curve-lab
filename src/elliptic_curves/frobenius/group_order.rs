use crate::elliptic_curves::frobenius::{CharacterSumPointCount, FrobeniusTrace, HasseInterval};
use crate::elliptic_curves::{CurveError, ExponentLowerBoundGroupOrderVerification};
use num_bigint::BigUint;

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
pub enum HasseGroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
}

impl HasseGroupOrderStrategy {
    pub fn as_group_order_strategy(self) -> GroupOrderStrategy {
        match self {
            Self::Auto => GroupOrderStrategy::Auto,
            Self::Exhaustive => GroupOrderStrategy::Exhaustive,
            Self::QuadraticCharacter => GroupOrderStrategy::QuadraticCharacter,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
    /// Tries to recover `#E(F_q)` from one known lower bound for the exponent
    /// `λ(E(F_q))` together with a Hasse-interval-producing base strategy.
    ///
    /// This route does not recount points from scratch. Instead it first runs
    /// the chosen `hasse_strategy` to obtain one Hasse interval `H(q)`, then
    /// asks whether that interval contains exactly one multiple of the
    /// supplied lower bound. It succeeds only in that uniqueness case.
    FromExponentLowerBoundAndPointCount {
        exponent_lower_bound: BigUint,
        hasse_strategy: HasseGroupOrderStrategy,
    },
}

/// Shared group-order result returned by curve-side methods.
///
/// This enum keeps the algorithmic route explicit instead of collapsing every
/// count immediately into a bare integer. Callers can still recover the common
/// arithmetic data such as `#E(F_q)`, the trace `t`, and the Hasse interval
/// without losing which strategy produced it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupOrderReport {
    ExhaustiveTrace(FrobeniusTrace),
    QuadraticCharacter(CharacterSumPointCount),
    FromExponentLowerBound(Box<ExponentLowerBoundGroupOrderVerification>),
}

impl GroupOrderReport {
    /// Returns the strategy used to build this report.
    pub fn strategy(&self) -> GroupOrderStrategy {
        match self {
            Self::ExhaustiveTrace(_) => GroupOrderStrategy::Exhaustive,
            Self::QuadraticCharacter(_) => GroupOrderStrategy::QuadraticCharacter,
            Self::FromExponentLowerBound(report) => {
                GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
                    exponent_lower_bound: report.exponent_lower_bound().clone(),
                    hasse_strategy: match report.group_order_report().strategy() {
                        GroupOrderStrategy::Auto => HasseGroupOrderStrategy::Auto,
                        GroupOrderStrategy::Exhaustive => HasseGroupOrderStrategy::Exhaustive,
                        GroupOrderStrategy::QuadraticCharacter => {
                            HasseGroupOrderStrategy::QuadraticCharacter
                        }
                        GroupOrderStrategy::FromExponentLowerBoundAndPointCount { .. } => {
                            unreachable!(
                                "the lower-bound report should wrap one base group-order strategy"
                            )
                        }
                    },
                }
            }
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => trace.field_order(),
            Self::QuadraticCharacter(report) => report.field_order(),
            Self::FromExponentLowerBound(report) => report.group_order_report().field_order(),
        }
    }

    /// Returns the counted curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => u128::from(trace.curve_order()),
            Self::QuadraticCharacter(report) => report.curve_order(),
            Self::FromExponentLowerBound(report) => report.group_order_report().curve_order(),
        }
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> i128 {
        match self {
            Self::ExhaustiveTrace(trace) => i128::from(trace.trace()),
            Self::QuadraticCharacter(report) => report.trace(),
            Self::FromExponentLowerBound(report) => report.group_order_report().trace(),
        }
    }

    /// Returns the discrete Hasse interval attached to the same `F_q`.
    pub fn hasse_interval(&self) -> HasseInterval {
        match self {
            Self::ExhaustiveTrace(trace) => trace.hasse_interval(),
            Self::QuadraticCharacter(report) => report.hasse_interval(),
            Self::FromExponentLowerBound(report) => report.group_order_report().hasse_interval(),
        }
    }

    /// Converts this report into the shared Frobenius-trace package.
    pub fn to_frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        match self {
            Self::ExhaustiveTrace(trace) => Ok(trace.clone()),
            Self::QuadraticCharacter(report) => report.to_frobenius_trace(),
            Self::FromExponentLowerBound(report) => {
                report.group_order_report().to_frobenius_trace()
            }
        }
    }
}
