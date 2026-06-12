use core::fmt;
use num_bigint::BigUint;

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
    /// A lower-bound-plus-Hasse route did not force one unique compatible
    /// group order.
    UnverifiedGroupOrderFromExponentLowerBound { lower_bound: BigUint },
    /// The supplied coefficients define a singular cubic.
    SingularCurve,
    /// The supplied affine coordinates do not satisfy the curve equation.
    PointNotOnCurve,
    /// Two function-field values were combined even though they belong to
    /// different short-Weierstrass curves.
    IncompatibleFunctionFieldCurves,
    /// A function-field element has zero norm and is therefore not invertible.
    NonInvertibleFunctionFieldElement,
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
            Self::UnverifiedGroupOrderFromExponentLowerBound { lower_bound } => {
                write!(
                    f,
                    "the chosen Hasse-based route does not verify one unique group order from the exponent lower bound {lower_bound}"
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
            Self::GroupAxiomViolation { axiom } => {
                write!(f, "finite-group axiom validation failed: {axiom}")
            }
        }
    }
}

impl std::error::Error for CurveError {}
