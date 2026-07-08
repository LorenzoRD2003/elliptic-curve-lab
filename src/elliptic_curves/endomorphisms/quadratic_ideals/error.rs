use core::fmt;

use crate::numerics::PositivePrimeError;

/// Failure modes for local prime behavior in an imaginary quadratic order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum QuadraticPrimeBehaviorError {
    InvalidPrime(PositivePrimeError),
    UnsupportedPrimeTwo,
}

impl fmt::Display for QuadraticPrimeBehaviorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPrime(error) => write!(f, "{error}"),
            Self::UnsupportedPrimeTwo => write!(
                f,
                "the first quadratic-prime-behavior API supports ℓ = 2 only when ℓ divides the conductor"
            ),
        }
    }
}

impl std::error::Error for QuadraticPrimeBehaviorError {}

impl From<PositivePrimeError> for QuadraticPrimeBehaviorError {
    fn from(error: PositivePrimeError) -> Self {
        Self::InvalidPrime(error)
    }
}
