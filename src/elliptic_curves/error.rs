use core::fmt;
use num_bigint::BigUint;

use crate::fields::FieldError;
use crate::polynomials::PolynomialError;

/// Errors returned when validating elliptic-curve models.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CurveError {
    /// The short Weierstrass model requires characteristic different from 2 and 3.
    UnsupportedCharacteristic { characteristic: u64 },
    /// A Frobenius-trace helper received invalid or unusable finite base-field metadata.
    InvalidFrobeniusBaseField {
        characteristic: u64,
        extension_degree: u32,
    },
    /// A Frobenius helper received a characteristic polynomial from a different finite base field.
    IncompatibleFrobeniusBaseField {
        curve_characteristic: u64,
        curve_extension_degree: u32,
        polynomial_characteristic: u64,
        polynomial_extension_degree: u32,
    },
    /// A Frobenius extension-count comparison received a curve over a finite field
    /// that is not compatible with the stored base field of the trace.
    IncompatibleFrobeniusTraceBaseField {
        trace_characteristic: u64,
        trace_extension_degree: u32,
        curve_characteristic: u64,
        curve_extension_degree: u32,
    },
    /// An isogeny Frobenius helper received domain and codomain curves over different finite fields.
    IncompatibleFrobeniusIsogenyBaseFields {
        domain_characteristic: u64,
        domain_extension_degree: u32,
        codomain_characteristic: u64,
        codomain_extension_degree: u32,
    },
    /// A Frobenius helper needs `[q]P`, but the current scalar-multiplication surface only accepts `u64`.
    UnsupportedFrobeniusFieldOrder { field_order: u128 },
    /// An absolute-Frobenius orbit helper was asked to iterate `π_p^k` on a curve
    /// that is not fixed by that Frobenius power.
    AbsoluteFrobeniusDoesNotPreserveCurve { power: u32 },
    /// A Frobenius-trace helper received an invalid or impossible curve order.
    InvalidCurveOrder { order: u64 },
    /// A Frobenius-trace helper received a trace that is incompatible with the requested `F_q`.
    InvalidFrobeniusTrace { trace: i64 },
    /// A Hasse-interval helper received an invalid or unsupported finite field order `q`.
    InvalidHasseIntervalFieldOrder { field_order: u128 },
    /// A character-sum point-count helper was asked to use a finite field
    /// whose quadratic-character route is not supported by the current backend.
    UnsupportedCharacterSumPointCount {
        characteristic: u64,
        extension_degree: u32,
    },
    /// A torsion helper received an invalid order parameter.
    InvalidTorsionOrder { order: usize },
    /// The odd-prime Schoof step requires an odd prime `ℓ` different from the
    /// field characteristic.
    InvalidSchoofOddPrime {
        odd_prime: usize,
        characteristic: u64,
    },
    /// The current Schoof route encountered one polynomial-domain failure
    /// while constructing or reducing an odd division polynomial.
    SchoofPolynomialFailure { error: PolynomialError },
    /// The current Schoof route unexpectedly asked the division-polynomial
    /// layer for an unsupported or ill-formed object after validating its
    /// odd-prime input.
    SchoofUnexpectedDivisionPolynomialModel,
    /// The current Schoof route unexpectedly reached one helper that still
    /// requires exhaustive base-field enumeration.
    SchoofUnexpectedFieldEnumerationRequirement,
    /// The current Schoof route unexpectedly reached one helper that still
    /// requires a square-root backend.
    SchoofUnexpectedSquareRootRequirement,
    /// The current automatic Schoof route stopped at one odd prime before the
    /// trace class became fully resolved.
    SchoofBlockedOnOddPrime { odd_prime: usize },
    /// The current automatic Schoof route produced a CRT trace class that is
    /// still too coarse for Hasse's bound to determine one unique trace.
    SchoofAmbiguousTraceClass { candidate_count: u128 },
    /// The current Schoof route produced intermediate data incompatible with
    /// Hasse's theorem.
    SchoofInconsistentWithHasse,
    /// An order-from-multiple helper received `M = 0`, which does not certify
    /// a finite point order.
    InvalidPointOrderMultiple { multiple: BigUint },
    /// An order-from-multiple helper received a factorization surface that is
    /// not a valid prime-power factorization of the supplied multiple.
    InvalidPointOrderMultipleFactorization { multiple: BigUint },
    /// An order-from-multiple helper was given an `M` such that `[M]P != O`.
    PointOrderMultipleDoesNotAnnihilatePoint { multiple: BigUint },
    /// A Hasse-interval search did not find any annihilating multiple.
    NoAnnihilatingMultipleInHasseInterval { lower: u128, upper: u128 },
    /// A reduced short-Weierstrass quotient `F[x, y] / (y^2 - f(x), g(x))`
    /// was asked to use the zero polynomial as the extra modulus `g(x)`.
    ZeroReducedCurveQuotientModulus,
    /// The deterministic group-order entry point was asked to run a strategy
    /// that needs a sampler-aware API.
    GroupOrderStrategyRequiresSampler { strategy: &'static str },
    /// Mestre's theorem in the current implementation is restricted to prime
    /// fields `F_p`, not extension fields `F_{p^r}` with `r > 1`.
    MestreRequiresPrimeField { extension_degree: u32 },
    /// The current Mestre route follows the prime-field theorem from the
    /// lecture notes, which assumes `p > 229`.
    MestrePrimeTooSmall { characteristic: u64 },
    /// The supplied sampler stopped before the Mestre route certified a
    /// unique group order.
    MestreSamplerExhausted,
    /// The requested Mestre iteration cap was reached before either side had
    /// a unique multiple in the Hasse interval.
    MestreIterationCapReached { max_iterations: usize },
    /// No genuinely quadratic twist was found for the current prime-field
    /// curve among the represented base-field factors.
    MestreQuadraticTwistUnavailable,
    /// The supplied coefficients define a singular cubic.
    SingularCurve,
    /// The supplied affine coordinates do not satisfy the curve equation.
    PointNotOnCurve,
    /// Two function-field values were combined even though they belong to
    /// different short-Weierstrass curves.
    IncompatibleFunctionFieldCurves,
    /// A function-field element has zero norm and is therefore not invertible.
    NonInvertibleFunctionFieldElement,
    /// A field-side validation or metadata query failed while serving a
    /// curve-side operation.
    Field(FieldError),
    /// An exhaustively checked finite-group axiom failed.
    GroupAxiomViolation { axiom: &'static str },
}

impl fmt::Display for CurveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedCharacteristic { characteristic } => write!(
                f,
                "short Weierstrass form requires characteristic different from 2 and 3, got {characteristic}"
            ),
            Self::InvalidFrobeniusBaseField {
                characteristic,
                extension_degree,
            } => write!(
                f,
                "invalid or unsupported finite base field for Frobenius trace computations: characteristic {characteristic}, extension degree {extension_degree}"
            ),
            Self::IncompatibleFrobeniusBaseField {
                curve_characteristic,
                curve_extension_degree,
                polynomial_characteristic,
                polynomial_extension_degree,
            } => write!(
                f,
                "Frobenius characteristic polynomial over F_{}^{} does not match curve base field F_{}^{}",
                polynomial_characteristic,
                polynomial_extension_degree,
                curve_characteristic,
                curve_extension_degree
            ),
            Self::IncompatibleFrobeniusTraceBaseField {
                trace_characteristic,
                trace_extension_degree,
                curve_characteristic,
                curve_extension_degree,
            } => write!(
                f,
                "Frobenius trace over F_{}^{} does not match curve base field F_{}^{} for extension-count comparison",
                trace_characteristic,
                trace_extension_degree,
                curve_characteristic,
                curve_extension_degree
            ),
            Self::IncompatibleFrobeniusIsogenyBaseFields {
                domain_characteristic,
                domain_extension_degree,
                codomain_characteristic,
                codomain_extension_degree,
            } => write!(
                f,
                "isogeny domain over F_{}^{} does not match codomain over F_{}^{} for Frobenius comparison",
                domain_characteristic,
                domain_extension_degree,
                codomain_characteristic,
                codomain_extension_degree
            ),
            Self::UnsupportedFrobeniusFieldOrder { field_order } => write!(
                f,
                "unsupported Frobenius field order for the current scalar-multiplication surface: {field_order}"
            ),
            Self::AbsoluteFrobeniusDoesNotPreserveCurve { power } => write!(
                f,
                "absolute Frobenius π_p^{power} does not preserve the current curve model"
            ),
            Self::InvalidCurveOrder { order } => {
                write!(
                    f,
                    "invalid curve order for Frobenius trace computations: {order}"
                )
            }
            Self::InvalidFrobeniusTrace { trace } => {
                write!(
                    f,
                    "invalid Frobenius trace for the requested finite base field: {trace}"
                )
            }
            Self::InvalidHasseIntervalFieldOrder { field_order } => {
                write!(
                    f,
                    "invalid or unsupported finite field order for Hasse interval computations: {field_order}"
                )
            }
            Self::UnsupportedCharacterSumPointCount {
                characteristic,
                extension_degree,
            } => write!(
                f,
                "character-sum point counting is not supported over F_{}^{} by the current quadratic-character backend",
                characteristic, extension_degree
            ),
            Self::InvalidTorsionOrder { order } => {
                write!(f, "torsion order must be a positive integer, got {order}")
            }
            Self::InvalidSchoofOddPrime {
                odd_prime,
                characteristic,
            } => write!(
                f,
                "the odd-prime Schoof step requires an odd prime l different from the field characteristic {characteristic}, got {odd_prime}"
            ),
            Self::SchoofPolynomialFailure { error } => write!(
                f,
                "the current Schoof route encountered a polynomial-domain failure: {error}"
            ),
            Self::SchoofUnexpectedDivisionPolynomialModel => write!(
                f,
                "the current Schoof route encountered an unexpected division-polynomial modeling failure after validating its odd-prime input"
            ),
            Self::SchoofUnexpectedFieldEnumerationRequirement => write!(
                f,
                "the current Schoof route unexpectedly fell back to a helper that requires exhaustive field enumeration"
            ),
            Self::SchoofUnexpectedSquareRootRequirement => write!(
                f,
                "the current Schoof route unexpectedly fell back to a helper that requires a square-root backend"
            ),
            Self::SchoofBlockedOnOddPrime { odd_prime } => write!(
                f,
                "the current Schoof route stopped at odd prime {odd_prime} before resolving the group order"
            ),
            Self::SchoofAmbiguousTraceClass { candidate_count } => write!(
                f,
                "the current Schoof CRT class still leaves {candidate_count} Hasse-compatible trace candidates"
            ),
            Self::SchoofInconsistentWithHasse => write!(
                f,
                "the current Schoof CRT class is incompatible with Hasse's theorem"
            ),
            Self::InvalidPointOrderMultiple { multiple } => {
                write!(
                    f,
                    "point-order-from-multiple requires a positive annihilating multiple, got {multiple}"
                )
            }
            Self::InvalidPointOrderMultipleFactorization { multiple } => {
                write!(
                    f,
                    "supplied factorization is not a valid prime-power factorization of the multiple {multiple}"
                )
            }
            Self::PointOrderMultipleDoesNotAnnihilatePoint { multiple } => {
                write!(
                    f,
                    "the supplied multiple does not annihilate the point: [{multiple}]P is not the identity"
                )
            }
            Self::NoAnnihilatingMultipleInHasseInterval { lower, upper } => {
                write!(
                    f,
                    "no annihilating multiple was found inside the Hasse interval [{lower}, {upper}]"
                )
            }
            Self::ZeroReducedCurveQuotientModulus => {
                write!(
                    f,
                    "the reduced short-Weierstrass quotient requires a non-zero univariate modulus"
                )
            }
            Self::GroupOrderStrategyRequiresSampler { strategy } => {
                write!(
                    f,
                    "the group-order strategy {strategy} requires `group_order_by_with_sampler(...)`"
                )
            }
            Self::MestreRequiresPrimeField { extension_degree } => {
                write!(
                    f,
                    "the current Mestre group-order route is implemented only over prime fields F_p, not extension degree {extension_degree}"
                )
            }
            Self::MestrePrimeTooSmall { characteristic } => {
                write!(
                    f,
                    "the current Mestre group-order route follows the p > 229 theorem and does not apply to characteristic {characteristic}"
                )
            }
            Self::MestreSamplerExhausted => {
                write!(
                    f,
                    "the supplied sampler exhausted before Mestre certified a unique group order"
                )
            }
            Self::MestreIterationCapReached { max_iterations } => {
                write!(
                    f,
                    "the Mestre group-order route hit its iteration cap of {max_iterations} before certifying a unique group order"
                )
            }
            Self::MestreQuadraticTwistUnavailable => {
                write!(
                    f,
                    "no genuinely quadratic twist could be selected for the current curve over its represented prime field"
                )
            }
            Self::SingularCurve => {
                write!(f, "short Weierstrass coefficients define a singular curve")
            }
            Self::PointNotOnCurve => {
                write!(f, "affine coordinates do not satisfy the curve equation")
            }
            Self::IncompatibleFunctionFieldCurves => {
                write!(f, "function-field elements belong to different curves")
            }
            Self::NonInvertibleFunctionFieldElement => {
                write!(f, "function-field element is not invertible")
            }
            Self::Field(error) => write!(f, "field error: {error}"),
            Self::GroupAxiomViolation { axiom } => {
                write!(f, "finite-group axiom validation failed: {axiom}")
            }
        }
    }
}

impl std::error::Error for CurveError {}

impl From<FieldError> for CurveError {
    fn from(error: FieldError) -> Self {
        Self::Field(error)
    }
}
