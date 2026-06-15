use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::isomorphisms::CurveIsomorphismError,
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Returns the short-Weierstrass model obtained from the scaling
    /// parameter `u`.
    pub fn scaled_by(&self, u: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&u).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }
        let u2 = F::square(&u);
        let u4 = F::square(&u2);
        let u6 = F::mul(&u4, &u2);
        Self::new(F::mul(&u4, self.a()), F::mul(&u6, self.b())).map_err(Into::into)
    }

    /// Returns whether `other` is exactly the short-Weierstrass model obtained
    /// by scaling this curve with the supplied parameter `u`.
    pub fn isomorphic_via_scale(&self, other: &Self, u: &F::Elem) -> bool {
        match self.scaled_by(u.clone()) {
            Ok(scaled_curve) => {
                F::eq(scaled_curve.a(), other.a()) && F::eq(scaled_curve.b(), other.b())
            }
            Err(_) => false,
        }
    }

    /// Returns the quadratic twist determined by the non-zero factor `d`.
    pub fn quadratic_twist(&self, d: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&d).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }
        let d2 = F::square(&d);
        let d3 = F::mul(&d2, &d);
        Self::new(F::mul(&d2, self.a()), F::mul(&d3, self.b())).map_err(Into::into)
    }
}
