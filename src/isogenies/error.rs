use core::fmt;

use crate::elliptic_curves::{CurveError, CurveIsomorphismError};

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
    /// Scalar multiplication by zero is not treated as an isogeny in the
    /// current educational surface.
    ///
    /// The map `[0]` is constant, so this crate keeps it out of the explicit
    /// isogeny constructors instead of silently modeling it alongside the
    /// non-constant multiplication-by-`n` maps.
    ZeroScalarIsNotIsogeny,
    /// An evaluated image point does not lie on the declared codomain curve.
    ///
    /// Exhaustive small-curve verifiers use this to report that some domain
    /// point was sent outside `E'(F_q)`, so the map does not even land on the
    /// claimed codomain.
    ImagePointNotOnCodomain,
    /// A declared kernel point failed to map to the codomain identity.
    ///
    /// This means the explicit kernel listing is not even contained in the
    /// actual kernel of the map.
    KernelPointDoesNotMapToIdentity,
    /// The map failed the additive homomorphism law on enumerated points.
    ///
    /// In other words, there exist points `P, Q` such that
    /// `phi(P + Q) != phi(P) + phi(Q)`.
    HomomorphismViolation,
    /// The explicit kernel listing does not coincide with the full fiber above
    /// the codomain identity.
    ///
    /// Exhaustive small-curve checks compute
    /// `{ P in E(F_q) : phi(P) = O }`
    /// directly and compare it against `kernel_points()`.
    KernelMismatch,
    /// Two explicit maps cannot be compared pointwise because they do not
    /// share the same concrete domain and codomain curves.
    ///
    /// This is distinct from a pointwise `false` result: before asking whether
    /// two maps agree on every domain point, both maps must first live between
    /// the same source and target curves.
    MapComparisonDomainCodomainMismatch,
    /// Two isogenies cannot be composed because the first codomain does not
    /// match the second domain.
    ///
    /// The upcoming composition layer will use this when attempting to
    /// form `psi ∘ phi` without a compatible middle curve.
    CompositionDomainCodomainMismatch,
    /// An exhaustive small-curve search failed to find a candidate dual
    /// isogeny.
    DualNotFound,
    /// A candidate dual was found, but the expected duality relations failed.
    ///
    /// In the small finite educational setting this means one of the checked
    /// identities, such as `hat(phi) ∘ phi = [n]`, did not hold.
    DualRelationViolation,
    /// The observed degree does not match the mathematically expected degree.
    DegreeMismatch,
    /// Two curves that should agree up to isomorphism did not admit a
    /// compatible base-field isomorphism in the attempted check.
    CurvesNotIsomorphic,
    /// A curve-isomorphism step failed in a way that does not map cleanly onto
    /// one of the existing isogeny-specific error categories.
    Isomorphism(CurveIsomorphismError),
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
            Self::ZeroScalarIsNotIsogeny => write!(
                formatter,
                "scalar multiplication by zero is not treated as an isogeny"
            ),
            Self::ImagePointNotOnCodomain => write!(
                formatter,
                "the isogeny sends an enumerated domain point outside the declared codomain"
            ),
            Self::KernelPointDoesNotMapToIdentity => write!(
                formatter,
                "a declared kernel point does not map to the codomain identity"
            ),
            Self::HomomorphismViolation => write!(
                formatter,
                "the isogeny violates the additive homomorphism law on enumerated points"
            ),
            Self::KernelMismatch => write!(
                formatter,
                "the explicit kernel points do not match the full identity fiber of the isogeny"
            ),
            Self::MapComparisonDomainCodomainMismatch => write!(
                formatter,
                "the two maps do not share the same concrete domain and codomain curves"
            ),
            Self::CompositionDomainCodomainMismatch => write!(
                formatter,
                "the isogeny codomain does not match the next isogeny domain for composition"
            ),
            Self::DualNotFound => write!(
                formatter,
                "no dual isogeny was found in the current exhaustive search"
            ),
            Self::DualRelationViolation => write!(
                formatter,
                "the candidate dual isogeny does not satisfy the expected duality relations"
            ),
            Self::DegreeMismatch => write!(
                formatter,
                "the isogeny degree does not match the expected degree"
            ),
            Self::CurvesNotIsomorphic => write!(
                formatter,
                "the compared curves are not isomorphic in the attempted check"
            ),
            Self::Isomorphism(error) => write!(
                formatter,
                "curve-isomorphism handling failed while building or comparing an isogeny: {error}"
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

impl From<CurveIsomorphismError> for IsogenyError {
    fn from(error: CurveIsomorphismError) -> Self {
        match error {
            CurveIsomorphismError::PointNotOnDomain => Self::Curve(CurveError::PointNotOnCurve),
            CurveIsomorphismError::ImagePointNotOnCodomain => Self::ImagePointNotOnCodomain,
            CurveIsomorphismError::CurvesNotIsomorphic => Self::CurvesNotIsomorphic,
            CurveIsomorphismError::UnsupportedCharacteristic { characteristic } => {
                Self::UnsupportedCharacteristic { characteristic }
            }
            CurveIsomorphismError::Curve(curve_error) => Self::Curve(curve_error),
            other => Self::Isomorphism(other),
        }
    }
}

impl std::error::Error for IsogenyError {}
