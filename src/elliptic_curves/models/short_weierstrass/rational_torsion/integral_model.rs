use num_bigint::BigInt;
use num_rational::BigRational;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::{
        isomorphisms::ShortWeierstrassIsomorphism, rational_torsion::RationalTorsionError,
    },
    traits::CurveIsomorphism,
};
use crate::fields::Q;
use crate::numerics::{lcm_biguint, rational_denominator_power_clearance};

/// Integral short-Weierstrass model prepared for rational-torsion search.
///
/// If the source curve is `E: y² = x³ + ax + b`, this witness wraps the
/// existing scaling isomorphism `ϕᵤ : E → E_int`, with
/// `ϕᵤ(x,y) = (u²x,u³y)`. The companion is the codomain
/// `E_int: y² = x³ + (u⁴a)x + u⁶b`, whose coefficients are integral.
#[derive(Clone, Debug)]
pub(crate) struct RationalIntegralModel {
    scaling: ShortWeierstrassIsomorphism<Q>,
}

impl RationalIntegralModel {
    /// Builds an integral model by choosing a positive integer scale.
    ///
    /// For each denominator `d`, the local scale contributes `p^⌈e/k⌉` when
    /// `p^e || d`, with `k = 4` for `a` and `k = 6` for `b`. This is the
    /// smallest positive integer forced by the two short-Weierstrass
    /// coefficient identities `a' = u⁴a` and `b' = u⁶b`.
    pub(crate) fn from_curve(
        source_curve: ShortWeierstrassCurve<Q>,
    ) -> Result<Self, RationalTorsionError> {
        let scale = integral_scale_for(&source_curve);
        Self::new(source_curve, scale)
    }

    pub(crate) fn new(
        source_curve: ShortWeierstrassCurve<Q>,
        scale: BigRational,
    ) -> Result<Self, RationalTorsionError> {
        let model = Self {
            scaling: ShortWeierstrassIsomorphism::new(source_curve, scale)?,
        };

        if !model.has_integral_coefficients() {
            return Err(RationalTorsionError::IntegralModelUnavailable);
        }

        Ok(model)
    }

    pub(crate) fn source_curve(&self) -> &ShortWeierstrassCurve<Q> {
        self.scaling.domain()
    }

    /// Returns the scaled companion model used by the rational-torsion search.
    pub(crate) fn curve(&self) -> &ShortWeierstrassCurve<Q> {
        self.scaling.codomain()
    }

    pub(crate) fn scale(&self) -> &BigRational {
        self.scaling.scaling_factor()
    }

    pub(crate) fn has_integral_coefficients(&self) -> bool {
        self.curve().a().is_integer() && self.curve().b().is_integer()
    }

    #[cfg(test)]
    pub(crate) fn to_integral_point(
        &self,
        point: &AffinePoint<Q>,
    ) -> Result<AffinePoint<Q>, RationalTorsionError> {
        self.scaling.evaluate(point).map_err(Into::into)
    }

    pub(crate) fn to_source_point(
        &self,
        point: &AffinePoint<Q>,
    ) -> Result<AffinePoint<Q>, RationalTorsionError> {
        self.scaling.inverse()?.evaluate(point).map_err(Into::into)
    }
}

fn integral_scale_for(curve: &ShortWeierstrassCurve<Q>) -> BigRational {
    let a_scale = rational_denominator_power_clearance(curve.a(), 4);
    let b_scale = rational_denominator_power_clearance(curve.b(), 6);
    BigRational::from_integer(BigInt::from(lcm_biguint(&a_scale, &b_scale)))
}

pub(super) fn integral_rational_to_bigint(
    value: &BigRational,
) -> Result<BigInt, RationalTorsionError> {
    if !value.is_integer() {
        return Err(RationalTorsionError::IntegralModelUnavailable);
    }

    Ok(value.to_integer())
}
