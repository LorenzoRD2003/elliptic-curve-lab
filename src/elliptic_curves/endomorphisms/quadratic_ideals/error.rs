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

/// Failure modes for constructing the first prime-norm ideal objects.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PrimeNormIdealError {
    InvalidPrime(PositivePrimeError),
    UnsupportedPrimeTwo,
    NonSplitPrime,
    NonInvertibleBecauseDividesConductor,
    RootDoesNotMatchPrimeBehavior,
}

impl fmt::Display for PrimeNormIdealError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPrime(error) => write!(f, "{error}"),
            Self::UnsupportedPrimeTwo => write!(
                f,
                "prime-norm ideal construction currently supports ℓ = 2 only through later dyadic ideal data"
            ),
            Self::NonSplitPrime => write!(
                f,
                "the first prime-norm ideal constructor only models split primes"
            ),
            Self::NonInvertibleBecauseDividesConductor => write!(
                f,
                "prime divides the conductor, so it is not invertible in this non-maximal order"
            ),
            Self::RootDoesNotMatchPrimeBehavior => write!(
                f,
                "chosen root is not one of the two split roots of Δ modulo ℓ"
            ),
        }
    }
}

impl std::error::Error for PrimeNormIdealError {}

impl From<QuadraticPrimeBehaviorError> for PrimeNormIdealError {
    fn from(error: QuadraticPrimeBehaviorError) -> Self {
        match error {
            QuadraticPrimeBehaviorError::InvalidPrime(error) => Self::InvalidPrime(error),
            QuadraticPrimeBehaviorError::UnsupportedPrimeTwo => Self::UnsupportedPrimeTwo,
        }
    }
}
