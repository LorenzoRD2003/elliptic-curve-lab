use crate::elliptic_curves::models::montgomery::MontgomeryCurve;
use crate::fields::traits::Field;

impl<F: Field> MontgomeryCurve<F> {
    /// Returns the discriminant `Δ = 16(A^2 - 4) / B^6`.
    pub fn discriminant(&self) -> F::Elem {
        let numerator = F::mul(
            &F::from_i64(16),
            &F::sub(&F::square(self.a()), &F::from_i64(4)),
        );
        let denominator = F::mul(&F::square(self.b()), &F::square(&F::square(self.b())));

        F::div(&numerator, &denominator)
            .expect("Montgomery discriminant denominator vanishes only when B = 0")
    }

    /// Returns the invariant `c4 = 16(A^2 - 3) / B^2`.
    pub fn c4(&self) -> F::Elem {
        let numerator = F::mul(
            &F::from_i64(16),
            &F::sub(&F::square(self.a()), &F::from_i64(3)),
        );
        let denominator = F::square(self.b());

        F::div(&numerator, &denominator)
            .expect("Montgomery c4 denominator vanishes only when B = 0")
    }

    /// Returns the invariant `c6 = 32A(9 - 2A^2) / B^3`.
    pub fn c6(&self) -> F::Elem {
        let numerator = F::mul(
            &F::from_i64(32),
            &F::mul(
                self.a(),
                &F::sub(
                    &F::from_i64(9),
                    &F::mul(&F::from_i64(2), &F::square(self.a())),
                ),
            ),
        );
        let denominator = F::mul(self.b(), &F::square(self.b()));

        F::div(&numerator, &denominator)
            .expect("Montgomery c6 denominator vanishes only when B = 0")
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    pub fn j_invariant(&self) -> F::Elem {
        let c4_cubed = F::cube(&self.c4());
        F::div(&c4_cubed, &self.discriminant())
            .expect("validated Montgomery curve has non-zero discriminant")
    }

    /// Returns whether this curve and `other` have the same `j`-invariant.
    pub fn has_same_j_invariant(&self, other: &Self) -> bool {
        F::eq(&self.j_invariant(), &other.j_invariant())
    }

    /// Returns the cubic right-hand side `x^3 + A x^2 + x`.
    pub(crate) fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        F::add(&F::add(&F::cube(x), &F::mul(self.a(), &F::square(x))), x)
    }
}
