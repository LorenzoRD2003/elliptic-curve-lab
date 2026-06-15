use core::fmt;

/// Failure modes for the canonical factorization `Δ = v^2 D_K`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum QuadraticDiscriminantFactorizationError {
    ZeroDiscriminant,
    PositiveDiscriminant,
    InvalidQuadraticOrderDiscriminant,
}

impl fmt::Display for QuadraticDiscriminantFactorizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDiscriminant => {
                write!(
                    f,
                    "quadratic discriminant factorization is undefined for Δ = 0"
                )
            }
            Self::PositiveDiscriminant => write!(
                f,
                "quadratic discriminant factorization currently expects an imaginary discriminant Δ < 0"
            ),
            Self::InvalidQuadraticOrderDiscriminant => write!(
                f,
                "quadratic discriminant is not congruent to 0 or 1 modulo 4, so it does not define a quadratic order"
            ),
        }
    }
}

impl std::error::Error for QuadraticDiscriminantFactorizationError {}

/// Failure modes for constructing an imaginary quadratic order
/// `O_f = Z + f O_K`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ImaginaryQuadraticOrderError {
    NonNegativeFundamentalDiscriminant,
    NonFundamentalDiscriminant,
    ZeroConductor,
    NonImaginaryOrderDiscriminant,
}

impl fmt::Display for ImaginaryQuadraticOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonNegativeFundamentalDiscriminant => write!(
                f,
                "imaginary quadratic orders require a negative fundamental discriminant D_K < 0"
            ),
            Self::NonFundamentalDiscriminant => write!(
                f,
                "imaginary quadratic order construction requires a fundamental discriminant D_K"
            ),
            Self::ZeroConductor => {
                write!(
                    f,
                    "imaginary quadratic order construction requires a positive conductor f >= 1"
                )
            }
            Self::NonImaginaryOrderDiscriminant => write!(
                f,
                "quadratic discriminant does not define an imaginary quadratic order"
            ),
        }
    }
}

impl std::error::Error for ImaginaryQuadraticOrderError {}

/// Failure modes for relative indices between imaginary quadratic orders.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum QuadraticOrderIndexError {
    DifferentQuadraticFields,
    NotSuborder,
}

impl fmt::Display for QuadraticOrderIndexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DifferentQuadraticFields => write!(
                f,
                "quadratic-order index is defined only for orders in the same imaginary quadratic field"
            ),
            Self::NotSuborder => write!(f, "quadratic-order index requires an inclusion O_f ⊆ O_g"),
        }
    }
}

impl std::error::Error for QuadraticOrderIndexError {}
