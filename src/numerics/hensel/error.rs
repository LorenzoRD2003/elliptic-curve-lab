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
    /// A general integer modulus must be at least `2`.
    TrivialModulus,
    /// The supplied root does not solve the required congruence modulo `p^k`.
    RootDoesNotSolveCurrentModulus,
    /// The derivative is not a unit modulo `p`, so the simple-root formula does
    /// not apply.
    SingularDerivativeModPrime,
    /// The odd-prime square-root helper received `p = 2`.
    EvenPrimeUnsupported,
    /// The current square-root modulo prime-power helper assumes `p` does not
    /// divide the radicand.
    RadicandDivisibleByPrimeUnsupported,
    /// A helper for radicands divisible by `p` received a unit radicand.
    RadicandNotDivisibleByPrime,
    /// The radicand is not a quadratic residue modulo the supplied prime.
    QuadraticNonResidueModPrime,
    /// The radicand has no square root modulo the requested prime power.
    NoSquareRootModuloPrimePower,
    /// The lifted residue did not certify an integer root inside the requested
    /// absolute bound.
    IntegerRootNotCertifiedInBound,
    /// The requested root bound would require a precision level larger than the
    /// current trace representation can store.
    TargetLevelOverflow,
    /// The configured simple-root seed scan would require enumerating too many
    /// residues modulo `p`.
    SeedScanLimitExceeded,
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
            Self::TrivialModulus => write!(
                formatter,
                "square roots modulo an integer require modulus at least 2"
            ),
            Self::RootDoesNotSolveCurrentModulus => write!(
                formatter,
                "the supplied root does not solve the expected congruence"
            ),
            Self::SingularDerivativeModPrime => write!(
                formatter,
                "simple-root Hensel lifting requires f'(x) to be non-zero modulo p"
            ),
            Self::EvenPrimeUnsupported => write!(
                formatter,
                "the odd-prime square-root helper does not handle p = 2"
            ),
            Self::RadicandDivisibleByPrimeUnsupported => write!(
                formatter,
                "square roots with radicand divisible by p are not supported by this helper yet"
            ),
            Self::RadicandNotDivisibleByPrime => write!(
                formatter,
                "this square-root helper expects the radicand to be divisible by p"
            ),
            Self::QuadraticNonResidueModPrime => write!(
                formatter,
                "the radicand is not a quadratic residue modulo p"
            ),
            Self::NoSquareRootModuloPrimePower => write!(
                formatter,
                "the radicand has no square root modulo the requested prime power"
            ),
            Self::IntegerRootNotCertifiedInBound => write!(
                formatter,
                "the lifted residue did not certify an integer root inside the requested bound"
            ),
            Self::TargetLevelOverflow => write!(
                formatter,
                "the requested Hensel precision level exceeds the supported range"
            ),
            Self::SeedScanLimitExceeded => write!(
                formatter,
                "the configured Hensel seed scan would enumerate too many residues"
            ),
        }
    }
}

impl std::error::Error for HenselLiftError {}
