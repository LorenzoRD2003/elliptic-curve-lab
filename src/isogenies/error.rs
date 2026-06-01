use core::fmt;

use crate::elliptic_curves::CurveError;

/// Errors produced while validating or constructing elliptic-curve isogenies.
///
/// The current isogeny work is intentionally educational, so this enum keeps
/// the intermediate mathematical failure modes explicit instead of collapsing
/// them into a single generic "invalid kernel" error. That makes it easier to
/// explain which subgroup axiom or curve-side hypothesis failed when building
/// a finite-kernel morphism.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyError {
    /// The proposed kernel has no points at all.
    ///
    /// A finite subgroup kernel must at least contain the identity, so an
    /// empty collection cannot define a valid finite-kernel isogeny.
    EmptyKernel,
    /// The proposed kernel does not contain the neutral element.
    ///
    /// For an additive elliptic-curve subgroup this means the set cannot be a
    /// subgroup, even before checking closure under inverses or addition.
    KernelDoesNotContainIdentity,
    /// At least one proposed kernel point is not a point of the domain curve.
    ///
    /// This captures the basic domain-compatibility requirement for kernels:
    /// every kernel point must live on the curve from which the isogeny is
    /// being constructed.
    KernelPointNotOnCurve,
    /// The proposed kernel is not stable under additive inverses.
    ///
    /// In other words, there exists a point `P` in the set whose inverse `-P`
    /// is missing, so the set is not a subgroup.
    KernelNotClosedUnderNegation,
    /// The proposed kernel is not closed under the elliptic-curve group law.
    ///
    /// This means there exist kernel points `P` and `Q` such that `P + Q` does
    /// not remain in the proposed kernel.
    KernelNotClosedUnderAddition,
    /// The current isogeny construction does not support the field
    /// characteristic.
    UnsupportedCharacteristic {
        /// Characteristic of the base field where the attempted construction
        /// lives.
        characteristic: u64,
    },
    /// A lower-level curve validation step failed while checking the domain,
    /// kernel points, or intermediate group operations.
    Curve(CurveError),
}

impl fmt::Display for IsogenyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyKernel => write!(
                formatter,
                "an isogeny kernel must contain at least the identity point"
            ),
            Self::KernelDoesNotContainIdentity => write!(
                formatter,
                "the proposed kernel does not contain the identity point"
            ),
            Self::KernelPointNotOnCurve => write!(
                formatter,
                "the proposed kernel contains a point that is not on the domain curve"
            ),
            Self::KernelNotClosedUnderNegation => write!(
                formatter,
                "the proposed kernel is not closed under point negation"
            ),
            Self::KernelNotClosedUnderAddition => write!(
                formatter,
                "the proposed kernel is not closed under the elliptic-curve group law"
            ),
            Self::UnsupportedCharacteristic { characteristic } => write!(
                formatter,
                "the current isogeny construction does not support characteristic {characteristic}"
            ),
            Self::Curve(error) => write!(
                formatter,
                "curve validation failed while building or evaluating an isogeny: {error}"
            ),
        }
    }
}

impl From<CurveError> for IsogenyError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl std::error::Error for IsogenyError {}
