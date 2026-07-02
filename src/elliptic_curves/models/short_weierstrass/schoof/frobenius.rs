use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::schoof::{ReducedCurveQuotient, ReducedEndomorphism},
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Builds the reduced Frobenius endomorphism
    /// `π_q = (x^q mod g(x), f(x)^((q-1)/2) mod g(x) * y)`.
    ///
    /// Here `E: y^2 = f(x)`, `q = #F_q`, and the reduction takes place
    /// in the quotient `F[x, y] / (y^2 - f(x), g(x))`.
    ///
    /// This is the odd-characteristic coordinate formula used in Schoof's
    /// algorithm after restricting to the reduced `x`-polynomial quotient.
    ///
    /// Complexity: if `m = deg g`, then `Θ((log q) m^2)` field operations,
    /// dominated by two modular exponentiations in `F[x] / (g(x))`.
    pub(crate) fn reduced_frobenius_endomorphism(
        &self,
        quotient: &ReducedCurveQuotient<F>,
    ) -> ReducedEndomorphism<F> {
        let x = DensePolynomial::new(vec![F::zero(), F::one()]);
        let cubic = self.to_cubic();
        let field_order = F::order().expect("finite field metadata should be valid");
        let x_map = DensePolynomial::pow_mod(&x, &field_order, quotient.modulus())
            .expect("the reduced-curve modulus is non-zero");
        let y_exponent = (&field_order - BigUint::from(1u8)) / BigUint::from(2u8);
        let y_scale = DensePolynomial::pow_mod(&cubic, &y_exponent, quotient.modulus())
            .expect("the reduced-curve modulus is non-zero");

        ReducedEndomorphism::new(quotient, x_map, y_scale)
    }
}

#[cfg(test)]
mod tests {

    use num_bigint::BigUint;

    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        short_weierstrass::schoof::{ReducedCurveQuotient, ReducedEndomorphism},
    };
    use crate::fields::traits::Field;
    use crate::polynomials::DensePolynomial;

    type F7 = crate::fields::Fp7;

    fn sample_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("sample curve should be smooth")
    }

    fn sample_quotient() -> ReducedCurveQuotient<F7> {
        ReducedCurveQuotient::new(
            sample_curve(),
            DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient")
    }

    #[test]
    fn reduced_frobenius_x_map_is_x_to_the_q_mod_g() {
        let curve = sample_curve();
        let quotient = sample_quotient();

        let frobenius = curve.reduced_frobenius_endomorphism(&quotient);
        let x = DensePolynomial::new(vec![F7::zero(), F7::one()]);
        let expected = DensePolynomial::pow_mod(&x, &BigUint::from(7u8), quotient.modulus())
            .expect("non-zero modulus should support modular exponentiation");

        assert_eq!(frobenius.x_map(), &expected);
    }

    #[test]
    fn reduced_frobenius_y_scale_is_f_to_the_q_minus_1_over_2_mod_g() {
        let curve = sample_curve();
        let quotient = sample_quotient();

        let frobenius = curve.reduced_frobenius_endomorphism(&quotient);
        let expected =
            DensePolynomial::pow_mod(&curve.to_cubic(), &BigUint::from(3u8), quotient.modulus())
                .expect("non-zero modulus should support modular exponentiation");

        assert_eq!(frobenius.y_scale(), &expected);
    }

    #[test]
    fn reduced_frobenius_is_reduced_modulo_g() {
        let curve = sample_curve();
        let quotient = sample_quotient();

        let frobenius: ReducedEndomorphism<F7> = curve.reduced_frobenius_endomorphism(&quotient);
        let modulus_degree = quotient
            .modulus()
            .degree()
            .expect("non-zero modulus has a degree");

        assert!(
            frobenius
                .x_map()
                .degree()
                .is_none_or(|degree| degree < modulus_degree)
        );
        assert!(
            frobenius
                .y_scale()
                .degree()
                .is_none_or(|degree| degree < modulus_degree)
        );
    }
}
use num_bigint::BigUint;
