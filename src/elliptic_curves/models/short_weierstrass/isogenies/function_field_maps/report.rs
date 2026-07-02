use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::{
        function_fields::ShortWeierstrassFunction,
        isogenies::function_field_maps::ShortWeierstrassFunctionFieldMap,
    },
};
use crate::fields::rational_function_field::RationalFunction;

/// Current separability classification surface for an explicit isogeny-like map.
///
/// The present implementation only certifies the separable case from the
/// differential multiplier. The remaining variants exist so the public API can
/// grow toward the full separable / purely inseparable / mixed story without a
/// breaking redesign.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IsogenySeparabilityKind {
    /// The differential multiplier is non-zero, so the map is certified separable.
    Separable,
    /// The map has been certified purely inseparable.
    ///
    /// This variant is reserved for later work.
    PurelyInseparable,
    /// The map has been certified as a non-trivial composition of separable and
    /// purely inseparable parts.
    ///
    /// This variant is reserved for later work.
    Mixed,
    /// The pullback data behaves like a constant map or otherwise does not
    /// model a non-constant elliptic-curve isogeny.
    ConstantOrInvalid,
    /// The current implementation cannot classify the zero-multiplier case more precisely.
    Unknown,
}

/// Differential pullback data for one short-Weierstrass function-field map.
///
/// Let
///
/// - `E : y^2 = x^3 + ax + b`
/// - `E' : y'^2 = x'^3 + a'x' + b'`
/// - `φ^* : F(E') -> F(E)`
///
/// be one stored pullback map with `X_φ = φ*(x')` and  `Y_φ = φ*(y')`.
///
/// The standard invariant differential on a short-Weierstrass curve is
///
/// `ω_E = dx / (2y)`.
///
/// In the present implementation, this report is intentionally modest:
///
/// - it fully certifies the separable case when `c_φ != 0`
/// - it records weaker classifications for the zero-multiplier case until the
///   broader inseparable-factor machinery is implemented
pub struct DifferentialPullbackReport<F: Field> {
    /// The source curve `E`, where the final comparison against `ω_E` lives.
    domain_curve: ShortWeierstrassCurve<F>,
    /// The target curve `E'`, whose invariant differential is being pulled back.
    codomain_curve: ShortWeierstrassCurve<F>,
    /// The function `X_φ = φ*(x') ∈ F(E)`.
    x_pullback: ShortWeierstrassFunction<F>,
    /// The function `Y_φ = φ*(y') ∈ F(E)`.
    y_pullback: ShortWeierstrassFunction<F>,
    /// The formal derivative `dX_φ/dx`, computed in the current function-field layer.
    dx_pullback: ShortWeierstrassFunction<F>,
    /// The factor multiplying `dx` in `φ*(ω_{E'}) = (dX_φ / (2Y_φ)) dx`.
    pulled_back_invariant_differential_factor: ShortWeierstrassFunction<F>,
    /// The function `c_φ = y (dX_φ/dx) / Y_φ` such that `φ*(ω_{E'}) = c_φ ω_E`.
    invariant_differential_multiplier: ShortWeierstrassFunction<F>,
    /// `Some(r)` when the computed multiplier visibly lies in `F(x) ⊂ F(E)`.
    rational_multiplier: Option<RationalFunction<F>>,
    /// The current separability-side classification inferred from the multiplier.
    separability_kind: IsogenySeparabilityKind,
}

impl<F: Field> DifferentialPullbackReport<F> {
    pub(crate) fn from_map_analysis(
        map: &ShortWeierstrassFunctionFieldMap<F>,
        dx_pullback: ShortWeierstrassFunction<F>,
        pulled_back_invariant_differential_factor: ShortWeierstrassFunction<F>,
        invariant_differential_multiplier: ShortWeierstrassFunction<F>,
        rational_multiplier: Option<RationalFunction<F>>,
        separability_kind: IsogenySeparabilityKind,
    ) -> Self
    where
        F::Elem: PartialEq,
    {
        Self {
            domain_curve: map.domain_curve().clone(),
            codomain_curve: map.codomain_curve().clone(),
            x_pullback: map.x_pullback().clone(),
            y_pullback: map.y_pullback().clone(),
            dx_pullback,
            pulled_back_invariant_differential_factor,
            invariant_differential_multiplier,
            rational_multiplier,
            separability_kind,
        }
    }

    pub fn domain_curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.domain_curve
    }

    pub fn codomain_curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.codomain_curve
    }

    pub fn x_pullback(&self) -> &ShortWeierstrassFunction<F> {
        &self.x_pullback
    }

    pub fn y_pullback(&self) -> &ShortWeierstrassFunction<F> {
        &self.y_pullback
    }

    pub fn dx_pullback(&self) -> &ShortWeierstrassFunction<F> {
        &self.dx_pullback
    }

    pub fn pulled_back_invariant_differential_factor(&self) -> &ShortWeierstrassFunction<F> {
        &self.pulled_back_invariant_differential_factor
    }

    pub fn invariant_differential_multiplier(&self) -> &ShortWeierstrassFunction<F> {
        &self.invariant_differential_multiplier
    }

    pub fn rational_multiplier(&self) -> Option<&RationalFunction<F>> {
        self.rational_multiplier.as_ref()
    }

    pub fn separability_kind(&self) -> IsogenySeparabilityKind {
        self.separability_kind
    }

    pub fn is_certified_separable(&self) -> bool {
        matches!(self.separability_kind, IsogenySeparabilityKind::Separable)
    }
}

impl<F: Field> fmt::Debug for DifferentialPullbackReport<F>
where
    F::Elem: fmt::Debug,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DifferentialPullbackReport")
            .field("domain_curve", &self.domain_curve)
            .field("codomain_curve", &self.codomain_curve)
            .field("x_pullback", &self.x_pullback)
            .field("y_pullback", &self.y_pullback)
            .field("dx_pullback", &self.dx_pullback)
            .field(
                "pulled_back_invariant_differential_factor",
                &self.pulled_back_invariant_differential_factor,
            )
            .field(
                "invariant_differential_multiplier",
                &self.invariant_differential_multiplier,
            )
            .field("rational_multiplier", &self.rational_multiplier)
            .field("separability_kind", &self.separability_kind)
            .finish()
    }
}
