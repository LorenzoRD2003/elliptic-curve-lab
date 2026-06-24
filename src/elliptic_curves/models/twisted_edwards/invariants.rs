use crate::elliptic_curves::TwistedEdwardsCurve;
use crate::fields::traits::Field;

impl<F: Field> TwistedEdwardsCurve<F> {
    /// /// Returns the Montgomery-normalized elliptic discriminant
    ///
    /// `Δ = a d (a - d)^4 / 16`.
    ///
    /// This normalization is the one induced by the canonical whole-curve
    /// bridge to the Montgomery model
    ///
    /// `A = 2(a + d)/(a - d)`, `B = 4/(a - d)`,
    ///
    /// together with the repo's Montgomery conventions.
    pub fn discriminant(&self) -> F::Elem {
        let a_minus_d = F::sub(self.a(), self.d());
        let numerator = F::mul(self.a(), &F::mul(self.d(), &F::pow(&a_minus_d, 4)));

        F::div(&numerator, &F::from_i64(16))
            .expect("characteristic different from 2 makes 16 invertible")
    }

    /// Returns the invariant  `c4 = a^2 + 14ad + d^2`.
    ///
    /// This is the Montgomery-normalized `c4` rewritten in the original
    /// twisted-Edwards coefficients.
    pub fn c4(&self) -> F::Elem {
        F::add(
            &F::add(
                &F::square(self.a()),
                &F::mul(&F::from_i64(14), &F::mul(self.a(), self.d())),
            ),
            &F::square(self.d()),
        )
    }

    /// Returns the invariant  `c6 = (a + d)(a^2 - 34ad + d^2)`.
    ///
    /// This is the Montgomery-normalized `c6` rewritten in the original
    /// twisted-Edwards coefficients.
    pub fn c6(&self) -> F::Elem {
        F::mul(
            &F::add(self.a(), self.d()),
            &F::add(
                &F::add(
                    &F::square(self.a()),
                    &F::mul(&F::from_i64(-34), &F::mul(self.a(), self.d())),
                ),
                &F::square(self.d()),
            ),
        )
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    ///
    /// i.e., `j = 16(a^2 + 14ad + d^2)^3 / (ad(a - d)^4)`
    pub fn j_invariant(&self) -> F::Elem {
        F::div(&F::cube(&self.c4()), &self.discriminant())
            .expect("validated twisted-Edwards curve has non-zero discriminant")
    }

    /// Returns whether this curve and `other` have the same `j`-invariant.
    pub fn has_same_j_invariant(&self, other: &Self) -> bool {
        F::eq(&self.j_invariant(), &other.j_invariant())
    }
}
