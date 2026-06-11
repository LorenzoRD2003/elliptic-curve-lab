use crate::elliptic_curves::isomorphisms::ShortWeierstrassIsomorphism;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::fields::EnumerableFiniteField;

impl<F: EnumerableFiniteField> ShortWeierstrassCurve<F> {
    /// Searches exhaustively for a base-field scaling isomorphism from this
    /// curve to `other`.
    pub fn find_isomorphism_to(&self, other: &Self) -> Option<ShortWeierstrassIsomorphism<F>> {
        for u in F::elements() {
            if F::is_zero(&u) {
                continue;
            }

            let scaled_curve = match self.scaled_by(u.clone()) {
                Ok(curve) => curve,
                Err(_) => continue,
            };

            if F::eq(scaled_curve.a(), other.a()) && F::eq(scaled_curve.b(), other.b()) {
                return ShortWeierstrassIsomorphism::new(
                    Self {
                        a: self.a.clone(),
                        b: self.b.clone(),
                    },
                    u,
                )
                .ok();
            }
        }

        None
    }

    /// Returns whether this curve is isomorphic to `other` over the current
    /// enumerable base field.
    pub fn is_isomorphic_to(&self, other: &Self) -> bool {
        self.find_isomorphism_to(other).is_some()
    }

    /// Returns every short-Weierstrass scaling automorphism of this curve over
    /// the current enumerable base field.
    pub fn automorphisms(&self) -> Vec<ShortWeierstrassIsomorphism<F>> {
        let mut automorphisms = Vec::new();
        for u in F::elements() {
            if !F::is_zero(&u)
                && let Ok(scaled_curve) = self.scaled_by(u.clone())
                && F::eq(scaled_curve.a(), self.a())
                && F::eq(scaled_curve.b(), self.b())
                && let Ok(isomorphism) = ShortWeierstrassIsomorphism::new(
                    Self {
                        a: self.a.clone(),
                        b: self.b.clone(),
                    },
                    u,
                )
            {
                automorphisms.push(isomorphism);
            }
        }
        automorphisms
    }
}
