use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::*;
use crate::polynomials::DensePolynomial;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Returns the cubic polynomial `x^3 + ax + b`.
    ///
    /// The coefficients are returned in ascending degree order, so the dense
    /// representation is `[b, a, 0, 1]`.
    pub fn to_cubic(&self) -> DensePolynomial<F> {
        DensePolynomial::<F>::new(vec![
            self.b().clone(),
            self.a().clone(),
            F::zero(),
            F::one(),
        ])
    }

    /// Returns the discriminant `Δ = -16(4a^3 + 27b^2)`.
    pub fn discriminant(&self) -> F::Elem {
        let four = F::from_i64(4);
        let minus_sixteen = F::from_i64(-16);
        let twenty_seven = F::from_i64(27);

        let four_a_cubed = F::mul(&four, &F::cube(self.a()));
        let twenty_seven_b_squared = F::mul(&twenty_seven, &F::square(self.b()));
        let inner = F::add(&four_a_cubed, &twenty_seven_b_squared);
        F::mul(&minus_sixteen, &inner)
    }

    /// Returns the classical Weierstrass invariant `c4 = -48a`.
    pub fn c4(&self) -> F::Elem {
        F::mul(&F::from_i64(-48), self.a())
    }

    /// Returns the classical Weierstrass invariant `c6 = -864b`.
    pub fn c6(&self) -> F::Elem {
        F::mul(&F::from_i64(-864), self.b())
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    pub fn j_invariant(&self) -> F::Elem {
        let c4_cubed = F::cube(&self.c4());
        F::div(&c4_cubed, &self.discriminant())
            .expect("validated short Weierstrass curve has non-zero discriminant")
    }

    /// Returns whether this curve and `other` have the same `j`-invariant.
    pub fn has_same_j_invariant(&self, other: &Self) -> bool {
        F::eq(&self.j_invariant(), &other.j_invariant())
    }

    /// Returns the cubic right-hand side `x^3 + ax + b`.
    pub(crate) fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        let x_cubed = F::cube(x);
        let ax = F::mul(self.a(), x);
        F::add(&F::add(&x_cubed, &ax), self.b())
    }
}
