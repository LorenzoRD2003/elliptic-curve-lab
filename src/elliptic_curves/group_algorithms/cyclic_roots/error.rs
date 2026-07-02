use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError, group_algorithms::cyclic_roots::input::CyclicPrimeRootInputError,
};

/// Errors returned by the staged cyclic prime-root algorithm.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CyclicPrimeRootError {
    /// The integer setup `(|G|, r)` is not a valid prime-root input.
    InvalidInput(CyclicPrimeRootInputError),
    /// A curve-side group operation failed.
    Curve(CurveError),
    /// The supplied `r`-Sylow generator `δ` does not have order `r^k`.
    InvalidSylowGenerator { expected_order: BigUint },
    /// The brute-force search did not find `x` with `α = xδ`.
    MissingSylowDiscreteLog { sylow_order: BigUint },
    /// The expected Bezout relation did not produce coprime data.
    MissingBezoutData {
        cofactor: BigUint,
        next_sylow_order: BigUint,
    },
}

impl From<CyclicPrimeRootInputError> for CyclicPrimeRootError {
    fn from(error: CyclicPrimeRootInputError) -> Self {
        Self::InvalidInput(error)
    }
}

impl From<CurveError> for CyclicPrimeRootError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}
