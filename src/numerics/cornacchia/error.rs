use core::fmt;

/// Failure modes for Cornacchia's algorithm with a supplied modular root.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CornacchiaError {
    /// The coefficient `d` in `x² + d y² = m` must be positive.
    ZeroCoefficient,
    /// The target integer `m` must be positive.
    ZeroTarget,
    /// The modular-root input must live modulo a non-trivial integer `m ≥ 2`.
    TrivialTarget,
    /// The supplied `r` is not a square root of `-d` modulo `m`.
    RootDoesNotSolveCongruence,
    /// Computing square roots of `-d` modulo `m` failed unexpectedly.
    ModularSquareRootFailure,
}

impl fmt::Display for CornacchiaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroCoefficient => write!(
                formatter,
                "Cornacchia's algorithm requires a positive coefficient d"
            ),
            Self::ZeroTarget => write!(
                formatter,
                "Cornacchia's algorithm requires a positive target m"
            ),
            Self::TrivialTarget => write!(
                formatter,
                "Cornacchia's algorithm with a supplied root requires m ≥ 2"
            ),
            Self::RootDoesNotSolveCongruence => {
                write!(formatter, "the supplied root must satisfy r² ≡ -d mod m")
            }
            Self::ModularSquareRootFailure => {
                write!(formatter, "computing square roots of -d modulo m failed")
            }
        }
    }
}

impl std::error::Error for CornacchiaError {}
