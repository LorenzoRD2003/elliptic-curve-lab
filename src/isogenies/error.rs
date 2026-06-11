use core::fmt;

use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::isomorphisms::CurveIsomorphismError;

/// Errors produced while validating the explicit finite kernel of an isogeny.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyKernelError {
    /// The proposed kernel has no points at all.
    EmptyKernel,
    /// The proposed kernel does not contain the neutral element.
    KernelDoesNotContainIdentity,
    /// At least one proposed kernel point is not a point of the domain curve.
    KernelPointNotOnCurve,
    /// The proposed kernel is not stable under additive inverses.
    KernelNotClosedUnderNegation,
    /// The proposed kernel is not closed under the elliptic-curve group law.
    KernelNotClosedUnderAddition,
}

impl fmt::Display for IsogenyKernelError {
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
        }
    }
}

/// Errors produced while constructing an isogeny that fails before kernel or
/// map verification even starts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyConstructionError {
    /// The current isogeny construction does not support the field
    /// characteristic.
    UnsupportedCharacteristic {
        /// Characteristic of the base field where the attempted construction
        /// lives.
        characteristic: u64,
    },
    /// Scalar multiplication by zero is not treated as an isogeny in the
    /// current educational surface.
    ZeroScalarIsNotIsogeny,
    /// Constructing Verschiebung from the direct `[p]^*` pullback requires
    /// inverting the absolute-Frobenius pullback on the corresponding
    /// coordinate.
    ///
    /// This records which coordinate failed that inversion.
    MissingInverseFrobeniusPreimageForVerschiebung {
        /// Which coordinate pullback failed to admit the needed inverse image.
        coordinate: &'static str,
    },
}

impl fmt::Display for IsogenyConstructionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedCharacteristic { characteristic } => write!(
                formatter,
                "the current isogeny construction does not support characteristic {characteristic}"
            ),
            Self::ZeroScalarIsNotIsogeny => write!(
                formatter,
                "scalar multiplication by zero is not treated as an isogeny"
            ),
            Self::MissingInverseFrobeniusPreimageForVerschiebung { coordinate } => write!(
                formatter,
                "the direct [p]^* pullback does not lie in the image of the absolute-Frobenius pullback for the {coordinate}-coordinate of Verschiebung"
            ),
        }
    }
}

/// Errors produced while verifying that an explicit map behaves like an
/// isogeny on small enumerated curves.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyVerificationError {
    /// An evaluated image point does not lie on the declared codomain curve.
    ImagePointNotOnCodomain,
    /// A declared kernel point failed to map to the codomain identity.
    KernelPointDoesNotMapToIdentity,
    /// The map failed the additive homomorphism law on enumerated points.
    HomomorphismViolation,
    /// The explicit kernel listing does not coincide with the full fiber above
    /// the codomain identity.
    KernelMismatch,
}

impl fmt::Display for IsogenyVerificationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
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
        }
    }
}

/// Errors produced while comparing, composing, or pulling back explicit
/// rational maps between short-Weierstrass function fields.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyMapError {
    /// Two explicit maps cannot be compared pointwise because they do not
    /// share the same concrete domain and codomain curves.
    MapComparisonDomainCodomainMismatch,
    /// Two isogenies cannot be composed because the first codomain does not
    /// match the second domain.
    CompositionDomainCodomainMismatch,
    /// The stored pullback functions do not live on the declared domain curve.
    FunctionFieldMapPullbackCurveMismatch,
    /// The stored pullback data does not satisfy the codomain equation after
    /// substitution.
    FunctionFieldMapCodomainEquationViolation,
    /// The function being pulled back does not live on the map's declared
    /// codomain curve.
    FunctionFieldMapSourceCurveMismatch,
    /// Substituting the stored `x`-pullback makes a rational denominator vanish
    /// identically in the domain function field.
    FunctionFieldMapDenominatorMapsToZero,
}

impl fmt::Display for IsogenyMapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MapComparisonDomainCodomainMismatch => write!(
                formatter,
                "the two maps do not share the same concrete domain and codomain curves"
            ),
            Self::CompositionDomainCodomainMismatch => write!(
                formatter,
                "the isogeny codomain does not match the next isogeny domain for composition"
            ),
            Self::FunctionFieldMapPullbackCurveMismatch => write!(
                formatter,
                "the stored function-field pullbacks do not live on the declared domain curve"
            ),
            Self::FunctionFieldMapCodomainEquationViolation => write!(
                formatter,
                "the stored function-field pullbacks do not satisfy the codomain curve equation"
            ),
            Self::FunctionFieldMapSourceCurveMismatch => write!(
                formatter,
                "the function being pulled back does not belong to the declared codomain function field"
            ),
            Self::FunctionFieldMapDenominatorMapsToZero => write!(
                formatter,
                "the pulled-back rational denominator becomes zero in the domain function field"
            ),
        }
    }
}

/// Errors produced while reasoning about dual isogenies or expected degrees.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DualIsogenyError {
    /// An exhaustive small-curve search failed to find a candidate dual
    /// isogeny.
    DualNotFound,
    /// A candidate dual was found, but the expected duality relations failed.
    DualRelationViolation,
    /// The observed degree does not match the mathematically expected degree.
    DegreeMismatch,
    /// Two curves that should agree up to isomorphism did not admit a
    /// compatible base-field isomorphism in the attempted check.
    CurvesNotIsomorphic,
}

impl fmt::Display for DualIsogenyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
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
        }
    }
}

/// Errors produced while checking a candidate Verschiebung against a chosen
/// Frobenius factorization.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VerschiebungError {
    /// A candidate Verschiebung pullback does not have source/target curves
    /// compatible with the chosen absolute Frobenius isogeny.
    DomainCodomainMismatch,
    /// The candidate Verschiebung failed the relation `V ∘ Frob_p = [p]`.
    LeftDualityViolation,
    /// The candidate Verschiebung failed the relation `Frob_p ∘ V = [p]`.
    RightDualityViolation,
}

impl fmt::Display for VerschiebungError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DomainCodomainMismatch => write!(
                formatter,
                "the candidate Verschiebung pullback does not match the source and target curves of the chosen Frobenius"
            ),
            Self::LeftDualityViolation => write!(
                formatter,
                "the candidate Verschiebung does not satisfy V ∘ Frob_p = [p]"
            ),
            Self::RightDualityViolation => write!(
                formatter,
                "the candidate Verschiebung does not satisfy Frob_p ∘ V = [p]"
            ),
        }
    }
}

/// Errors produced while validating or constructing elliptic-curve isogenies.
///
/// The current isogeny work is intentionally educational, so the top-level
/// error type keeps the main mathematical failure families explicit instead of
/// flattening everything into one large enum.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsogenyError {
    /// Explicit finite-kernel validation failed.
    Kernel(IsogenyKernelError),
    /// A construction policy or characteristic restriction failed.
    Construction(IsogenyConstructionError),
    /// Exhaustive small-curve verification failed.
    Verification(IsogenyVerificationError),
    /// Map comparison, composition, or function-field pullback data failed.
    Map(IsogenyMapError),
    /// Duality or degree reasoning failed.
    Dual(DualIsogenyError),
    /// A candidate Verschiebung failed one of its certified checks.
    Verschiebung(VerschiebungError),
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
            Self::Kernel(error) => write!(formatter, "{error}"),
            Self::Construction(error) => write!(formatter, "{error}"),
            Self::Verification(error) => write!(formatter, "{error}"),
            Self::Map(error) => write!(formatter, "{error}"),
            Self::Dual(error) => write!(formatter, "{error}"),
            Self::Verschiebung(error) => write!(formatter, "{error}"),
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

impl From<IsogenyKernelError> for IsogenyError {
    fn from(error: IsogenyKernelError) -> Self {
        Self::Kernel(error)
    }
}

impl From<IsogenyConstructionError> for IsogenyError {
    fn from(error: IsogenyConstructionError) -> Self {
        Self::Construction(error)
    }
}

impl From<IsogenyVerificationError> for IsogenyError {
    fn from(error: IsogenyVerificationError) -> Self {
        Self::Verification(error)
    }
}

impl From<IsogenyMapError> for IsogenyError {
    fn from(error: IsogenyMapError) -> Self {
        Self::Map(error)
    }
}

impl From<DualIsogenyError> for IsogenyError {
    fn from(error: DualIsogenyError) -> Self {
        Self::Dual(error)
    }
}

impl From<VerschiebungError> for IsogenyError {
    fn from(error: VerschiebungError) -> Self {
        Self::Verschiebung(error)
    }
}

impl From<CurveIsomorphismError> for IsogenyError {
    fn from(error: CurveIsomorphismError) -> Self {
        Self::Isomorphism(error)
    }
}

impl From<CurveError> for IsogenyError {
    fn from(error: CurveError) -> Self {
        Self::Curve(error)
    }
}

impl std::error::Error for IsogenyError {}
