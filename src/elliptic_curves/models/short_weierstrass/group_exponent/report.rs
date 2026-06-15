use num_bigint::BigUint;

use crate::elliptic_curves::short_weierstrass::group_exponent::{
    ExponentAccumulationReport, GroupExponentStrategy,
};

/// Shared group-exponent report returned by the unified curve-side API.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupExponentReport<P> {
    Exhaustive(BigUint),
    RandomPoints(Box<ExponentAccumulationReport<P>>),
}

impl<P> GroupExponentReport<P> {
    /// Returns the strategy used to build this report.
    pub fn strategy(&self) -> GroupExponentStrategy {
        match self {
            Self::Exhaustive(_) => GroupExponentStrategy::Exhaustive,
            Self::RandomPoints(report) => report.strategy(),
        }
    }

    /// Returns the best current lower bound for the exponent.
    ///
    /// For the exhaustive route this is exact. For the random-point route this
    /// is the running `lcm` lower bound accumulated from sampled point orders.
    pub fn exponent_lower_bound(&self) -> &BigUint {
        match self {
            Self::Exhaustive(exponent) => exponent,
            Self::RandomPoints(report) => report.exponent_lower_bound(),
        }
    }

    /// Returns the exact exponent when the chosen route computes it
    /// exhaustively.
    pub fn exact_exponent(&self) -> Option<&BigUint> {
        match self {
            Self::Exhaustive(exponent) => Some(exponent),
            Self::RandomPoints(_) => None,
        }
    }
}
