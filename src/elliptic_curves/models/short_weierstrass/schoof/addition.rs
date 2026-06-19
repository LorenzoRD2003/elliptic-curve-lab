use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::schoof::{
        QuotientInverseResult, ReducedCurveQuotient, ReducedEndomorphism,
        ReducedEndomorphismAdditiveResult,
    },
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Adds two reduced endomorphisms of the form `(a(x), b(x) y)`.
    ///
    /// When the stored `x`-images differ, this uses the secant slope
    ///
    /// `r = (b1 - b2) / (a1 - a2)`.
    ///
    /// When both reduced maps are literally equal, it uses the tangent slope
    ///
    /// `r = (3a^2 + A) / (2bf)`.
    ///
    /// When the stored `x`-images agree and the stored `y`-scales are
    /// additive inverses, this returns the additive zero endomorphism.
    ///
    /// If the required denominator is not a unit modulo `g(x)`, this returns
    /// a gcd witness through
    /// [`ReducedEndomorphismAdditiveResult::NonUnitDenominator`].
    ///
    /// Complexity: `Θ(m^2)` field operations, where `m = deg g`.
    pub(crate) fn add_reduced_endomorphisms(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        lhs: &ReducedEndomorphism<F>,
        rhs: &ReducedEndomorphism<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        if lhs == rhs {
            return self.double_reduced_endomorphism(quotient, lhs);
        }
        if self.reduced_endomorphisms_are_additive_inverses(lhs, rhs) {
            return ReducedEndomorphismAdditiveResult::Zero;
        }

        let numerator = lhs.y_scale().sub(rhs.y_scale());
        let denominator = lhs.x_map().sub(rhs.x_map());
        self.finish_reduced_endomorphism_sum_from_fraction(
            quotient,
            lhs,
            rhs,
            numerator,
            denominator,
        )
    }

    /// Doubles one reduced endomorphism of the form `(a(x), b(x) y)`.
    ///
    /// This uses the tangent slope
    ///
    /// `r = (3a^2 + A) / (2bf)`.
    ///
    /// If `b(x) = 0`, then the image lies on the `y = 0` branch and the
    /// tangent is vertical, so the double is the additive zero endomorphism.
    ///
    /// If the denominator is not a unit modulo `g(x)`, this returns a gcd
    /// witness through
    /// [`ReducedEndomorphismAdditiveResult::NonUnitDenominator`].
    ///
    /// Complexity: `Θ(m^2)` field operations, where `m = deg g`.
    pub(crate) fn double_reduced_endomorphism(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        value: &ReducedEndomorphism<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        if value.y_scale().is_zero() {
            return ReducedEndomorphismAdditiveResult::Zero;
        }

        let x_squared = value.x_map().mul(value.x_map());
        let three_x_squared = x_squared.scale(&F::from_i64(3));
        let numerator = three_x_squared.add(&DensePolynomial::constant(self.a().clone()));
        let denominator = value.y_scale().mul(&self.to_cubic()).scale(&F::from_i64(2));

        self.finish_reduced_endomorphism_sum_from_fraction(
            quotient,
            value,
            value,
            numerator,
            denominator,
        )
    }

    fn finish_reduced_endomorphism_sum_from_fraction(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        lhs: &ReducedEndomorphism<F>,
        rhs: &ReducedEndomorphism<F>,
        numerator: DensePolynomial<F>,
        denominator: DensePolynomial<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        match self.try_reduce_fraction_mod_g(quotient, numerator, denominator) {
            Ok(slope) => ReducedEndomorphismAdditiveResult::Value(
                self.reduced_endomorphism_sum_from_slope(quotient, lhs, rhs, slope),
            ),
            Err(witness_gcd) => {
                ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd }
            }
        }
    }

    fn reduced_endomorphisms_are_additive_inverses(
        &self,
        lhs: &ReducedEndomorphism<F>,
        rhs: &ReducedEndomorphism<F>,
    ) -> bool {
        lhs.x_map() == rhs.x_map() && lhs.y_scale() == &rhs.y_scale().neg()
    }

    fn try_reduce_fraction_mod_g(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        numerator: DensePolynomial<F>,
        denominator: DensePolynomial<F>,
    ) -> Result<DensePolynomial<F>, DensePolynomial<F>> {
        match quotient.try_invert_poly(&denominator) {
            QuotientInverseResult::Inverse(inverse) => {
                Ok(quotient.reduce_poly(&numerator.mul(&inverse)))
            }
            QuotientInverseResult::NonUnit { witness_gcd } => Err(witness_gcd),
        }
    }

    fn reduced_endomorphism_sum_from_slope(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        lhs: &ReducedEndomorphism<F>,
        rhs: &ReducedEndomorphism<F>,
        slope: DensePolynomial<F>,
    ) -> ReducedEndomorphism<F> {
        let cubic = self.to_cubic();
        let x_map = quotient.reduce_poly(
            &slope
                .mul(&slope)
                .mul(&cubic)
                .sub(lhs.x_map())
                .sub(rhs.x_map()),
        );
        let y_scale = quotient.reduce_poly(&slope.mul(&lhs.x_map().sub(&x_map)).sub(lhs.y_scale()));
        ReducedEndomorphism::new(quotient, x_map, y_scale)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        short_weierstrass::schoof::{
            ReducedCurveQuotient, ReducedEndomorphism, ReducedEndomorphismAdditiveResult,
        },
    };
    use crate::fields::{Fp, traits::Field};
    use crate::polynomials::DensePolynomial;

    type F7 = Fp<7>;

    fn sample_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("sample curve should be smooth")
    }

    fn sample_quotient() -> ReducedCurveQuotient<F7> {
        ReducedCurveQuotient::new(
            sample_curve(),
            DensePolynomial::new(vec![F7::from_i64(-1), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient")
    }

    fn invertible_doubling_quotient() -> ReducedCurveQuotient<F7> {
        ReducedCurveQuotient::new(
            sample_curve(),
            DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient")
    }

    #[test]
    fn successful_endomorphism_sum_returns_the_reduced_coordinate_formula() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let lhs = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::zero(), F7::one()]),
            DensePolynomial::constant(F7::one()),
        );
        let rhs = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::one(), F7::one()]),
            DensePolynomial::constant(F7::one()),
        );

        let ReducedEndomorphismAdditiveResult::Value(sum) =
            curve.add_reduced_endomorphisms(&quotient, &lhs, &rhs)
        else {
            panic!("the chosen denominator is a non-zero constant and should be invertible");
        };

        assert_eq!(
            sum.x_map(),
            &DensePolynomial::new(vec![F7::from_i64(-1), F7::from_i64(-2)])
        );
        assert_eq!(sum.y_scale(), &DensePolynomial::constant(F7::from_i64(-1)));
    }

    #[test]
    fn non_unit_secant_denominator_returns_its_gcd_witness() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let lhs = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::zero(), F7::one()]),
            DensePolynomial::constant(F7::one()),
        );
        let rhs = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::constant(F7::one()),
            DensePolynomial::constant(F7::from_i64(2)),
        );

        let ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd } =
            curve.add_reduced_endomorphisms(&quotient, &lhs, &rhs)
        else {
            panic!("x - 1 should not be invertible modulo x^2 - 1");
        };

        assert_eq!(
            witness_gcd,
            DensePolynomial::new(vec![F7::from_i64(-1), F7::one()])
        );
    }

    #[test]
    fn doubling_helper_matches_addition_of_equal_maps() {
        let curve = sample_curve();
        let quotient = invertible_doubling_quotient();
        let value = ReducedEndomorphism::identity(&quotient);

        let doubled = curve.double_reduced_endomorphism(&quotient, &value);
        let added = curve.add_reduced_endomorphisms(&quotient, &value, &value);

        match (doubled, added) {
            (
                ReducedEndomorphismAdditiveResult::Value(doubled),
                ReducedEndomorphismAdditiveResult::Value(added),
            ) => {
                assert_eq!(doubled.x_map(), added.x_map());
                assert_eq!(doubled.y_scale(), added.y_scale());
            }
            _ => panic!(
                "identity should double through an invertible tangent denominator in this quotient"
            ),
        }
    }

    #[test]
    fn opposite_reduced_endomorphisms_add_to_the_additive_zero() {
        let curve = sample_curve();
        let quotient = invertible_doubling_quotient();
        let value = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(6)]),
        );
        let inverse = value.additive_inverse();

        match curve.add_reduced_endomorphisms(&quotient, &value, &inverse) {
            ReducedEndomorphismAdditiveResult::Zero => {}
            other => panic!("expected additive zero for opposite maps, got {other:?}"),
        }
    }
}
