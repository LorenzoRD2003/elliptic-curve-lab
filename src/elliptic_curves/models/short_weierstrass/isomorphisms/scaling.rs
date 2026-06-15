use core::fmt;
use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::CurveIsomorphismError,
    traits::{CurveIsomorphism, CurveModel},
};
use crate::fields::traits::{EnumerableFiniteField, Field};

/// Explicit short-Weierstrass scaling isomorphism determined by a parameter `u`.
///
/// The intended mathematical convention is the map `ϕ_u : E -> E'` defined on
/// affine points by `(x, y) -> (u^2 x, u^3 y)`.
///
/// If the domain curve is written in short-Weierstrass form as
/// `E: y^2 = x^3 + ax + b`, then the image curve is `E': y^2 = x^3 + a'x + b'`
/// with transformed coefficients `a' = u^4 a`, `b' = u^6 b`.
///
/// This type treats `domain` and the parameter `u` as the primary data. The
/// codomain is derived automatically from those values instead of being stored
/// as a second source of truth.
#[derive(Clone)]
pub struct ShortWeierstrassIsomorphism<F: Field> {
    domain: ShortWeierstrassCurve<F>,
    codomain: ShortWeierstrassCurve<F>,
    u: F::Elem,
}

impl<F: Field> fmt::Debug for ShortWeierstrassIsomorphism<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShortWeierstrassIsomorphism")
            .field("domain", &self.domain)
            .field("codomain", &self.codomain)
            .field("u", &self.u)
            .finish()
    }
}

impl<F: Field> ShortWeierstrassIsomorphism<F> {
    /// Builds the short-Weierstrass scaling isomorphism determined by `u`.
    ///
    /// Since `ϕ_u` must be invertible, `u` must be invertible in the base field.
    pub fn new(
        domain: ShortWeierstrassCurve<F>,
        u: F::Elem,
    ) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&u).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }
        let codomain = Self::derive_codomain_from(&domain, &u)?;
        Ok(Self {
            domain,
            codomain,
            u,
        })
    }

    /// Returns the scaling parameter `u` from `ϕ_u`.
    pub fn scaling_factor(&self) -> &F::Elem {
        &self.u
    }

    /// Returns the inverse isomorphism `ϕ_u^{-1}`.
    ///
    /// Under the convention `ϕ_u(x, y) = (u^2 x, u^3 y)`, the inverse map
    /// is the short-Weierstrass scaling determined by `u^{-1}` from the
    /// derived codomain `E'` back to the original domain `E`.
    pub fn inverse(&self) -> Result<Self, CurveIsomorphismError> {
        let inverse_u = F::inv(&self.u).ok_or(CurveIsomorphismError::NonInvertibleScale)?;
        Self::new(self.codomain.clone(), inverse_u)
    }

    fn derive_codomain_from(
        domain: &ShortWeierstrassCurve<F>,
        u: &F::Elem,
    ) -> Result<ShortWeierstrassCurve<F>, CurveIsomorphismError> {
        let u2 = F::square(u);
        let u4 = F::square(&u2);
        let u6 = F::mul(&u4, &u2);
        ShortWeierstrassCurve::new(F::mul(&u4, domain.a()), F::mul(&u6, domain.b()))
            .map_err(Into::into)
    }
}

impl<F: Field + EnumerableFiniteField + Clone> ShortWeierstrassCurve<F>
where
    F::Elem: Clone + Eq + Hash,
{
    /// Enumerates all short-Weierstrass base-field scaling isomorphisms from
    /// `self` onto `target` by exhaustive search over the represented field.
    ///
    /// This is intentionally crate-private and educational: it is a tiny-field
    /// witness search used by current short-Weierstrass dual-isogeny routines.
    pub(crate) fn exhaustive_isomorphisms_to(
        &self,
        target: &ShortWeierstrassCurve<F>,
    ) -> Vec<ShortWeierstrassIsomorphism<F>> {
        F::elements()
            .into_iter()
            .filter_map(|u| ShortWeierstrassIsomorphism::new(self.clone(), u).ok())
            .filter(|isomorphism| isomorphism.codomain() == target)
            .collect()
    }
}

impl<F: Field> CurveIsomorphism for ShortWeierstrassIsomorphism<F> {
    type Domain = ShortWeierstrassCurve<F>;
    type Codomain = ShortWeierstrassCurve<F>;

    fn domain(&self) -> &Self::Domain {
        &self.domain
    }

    fn codomain(&self) -> &Self::Codomain {
        &self.codomain
    }

    fn evaluate(
        &self,
        point: &<Self::Domain as CurveModel>::Point,
    ) -> Result<<Self::Codomain as CurveModel>::Point, CurveIsomorphismError> {
        if !self.domain.contains(point) {
            return Err(CurveIsomorphismError::PointNotOnDomain);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::infinity()),
            AffinePoint::Finite { x, y } => {
                let image =
                    AffinePoint::new(F::mul(&F::square(&self.u), x), F::mul(&F::cube(&self.u), y));

                if !self.codomain.contains(&image) {
                    return Err(CurveIsomorphismError::ImagePointNotOnCodomain);
                }

                Ok(image)
            }
        }
    }
}
