use core::fmt;

/// Failure modes for the basic simple-root Hensel lifting helper.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum HenselLiftError {
    /// The target precision must be at least modulo `p`.
    ZeroTargetLevel,
    /// A polynomial must contain at least one coefficient.
    EmptyPolynomial,
    /// The simple-root lifting helper needs a non-constant polynomial.
    ConstantPolynomial,
    /// The modulus parameter must be a prime integer.
    NonPrimeModulus,
    /// The supplied root does not solve the required congruence modulo `p^k`.
    RootDoesNotSolveCurrentModulus,
    /// The derivative is not a unit modulo `p`, so the simple-root formula does
    /// not apply.
    SingularDerivativeModPrime,
}

impl fmt::Display for HenselLiftError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroTargetLevel => {
                write!(formatter, "Hensel lifting requires target level at least 1")
            }
            Self::EmptyPolynomial => write!(formatter, "Hensel lifting requires a polynomial"),
            Self::ConstantPolynomial => write!(
                formatter,
                "simple-root Hensel lifting requires a non-constant polynomial"
            ),
            Self::NonPrimeModulus => write!(
                formatter,
                "Hensel lifting currently requires a prime modulus"
            ),
            Self::RootDoesNotSolveCurrentModulus => write!(
                formatter,
                "the supplied root does not solve the expected congruence"
            ),
            Self::SingularDerivativeModPrime => write!(
                formatter,
                "simple-root Hensel lifting requires f'(x) to be non-zero modulo p"
            ),
        }
    }
}

impl std::error::Error for HenselLiftError {}
