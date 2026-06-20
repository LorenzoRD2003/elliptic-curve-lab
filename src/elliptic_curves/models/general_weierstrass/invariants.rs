use crate::elliptic_curves::models::general_weierstrass::GeneralWeierstrassCurve;
use crate::fields::traits::Field;

impl<F: Field> GeneralWeierstrassCurve<F> {
    /// Returns the invariant `b2 = a1^2 + 4a2`.
    pub fn b2(&self) -> F::Elem {
        F::add(&F::square(self.a1()), &F::mul(&F::from_i64(4), self.a2()))
    }

    /// Returns the invariant `b4 = 2a4 + a1a3`.
    pub fn b4(&self) -> F::Elem {
        F::add(
            &F::mul(&F::from_i64(2), self.a4()),
            &F::mul(self.a1(), self.a3()),
        )
    }

    /// Returns the invariant `b6 = a3^2 + 4a6`.
    pub fn b6(&self) -> F::Elem {
        F::add(&F::square(self.a3()), &F::mul(&F::from_i64(4), self.a6()))
    }

    /// Returns the invariant
    ///
    /// `b8 = a1^2 a6 + 4a2 a6 - a1 a3 a4 + a2 a3^2 - a4^2`.
    pub fn b8(&self) -> F::Elem {
        let a1_squared_a6 = F::mul(&F::square(self.a1()), self.a6());
        let four_a2_a6 = F::mul(&F::from_i64(4), &F::mul(self.a2(), self.a6()));
        let a1_a3_a4 = F::mul(&F::mul(self.a1(), self.a3()), self.a4());
        let a2_a3_squared = F::mul(self.a2(), &F::square(self.a3()));
        let a4_squared = F::square(self.a4());

        F::sub(
            &F::add(&F::add(&a1_squared_a6, &four_a2_a6), &a2_a3_squared),
            &F::add(&a1_a3_a4, &a4_squared),
        )
    }

    /// Returns the invariant `c4 = b2^2 - 24b4`.
    pub fn c4(&self) -> F::Elem {
        F::sub(
            &F::square(&self.b2()),
            &F::mul(&F::from_i64(24), &self.b4()),
        )
    }

    /// Returns the invariant `c6 = -b2^3 + 36b2b4 - 216b6`.
    pub fn c6(&self) -> F::Elem {
        let minus_b2_cubed = F::neg(&F::cube(&self.b2()));
        let thirty_six_b2_b4 = F::mul(&F::from_i64(36), &F::mul(&self.b2(), &self.b4()));
        let two_hundred_sixteen_b6 = F::mul(&F::from_i64(216), &self.b6());

        F::sub(
            &F::add(&minus_b2_cubed, &thirty_six_b2_b4),
            &two_hundred_sixteen_b6,
        )
    }

    /// Returns the discriminant
    ///
    /// `Δ = -b2^2 b8 - 8b4^3 - 27b6^2 + 9b2b4b6`.
    pub fn discriminant(&self) -> F::Elem {
        let minus_b2_squared_b8 = F::neg(&F::mul(&F::square(&self.b2()), &self.b8()));
        let minus_eight_b4_cubed = F::neg(&F::mul(&F::from_i64(8), &F::cube(&self.b4())));
        let minus_twenty_seven_b6_squared =
            F::neg(&F::mul(&F::from_i64(27), &F::square(&self.b6())));
        let nine_b2_b4_b6 = F::mul(
            &F::from_i64(9),
            &F::mul(&F::mul(&self.b2(), &self.b4()), &self.b6()),
        );

        F::add(
            &F::add(&minus_b2_squared_b8, &minus_eight_b4_cubed),
            &F::add(&minus_twenty_seven_b6_squared, &nine_b2_b4_b6),
        )
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    pub fn j_invariant(&self) -> F::Elem {
        let c4_cubed = F::cube(&self.c4());
        F::div(&c4_cubed, &self.discriminant())
            .expect("validated general Weierstrass curve has non-zero discriminant")
    }
}
