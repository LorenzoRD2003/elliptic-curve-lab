use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::function_fields::{ShortWeierstrassFunction, ShortWeierstrassFunctionField},
    short_weierstrass::schoof::{
        QuotientInverseResult, ReducedCurveQuotient, ReducedEndomorphism,
        ReducedEndomorphismAdditiveResult,
    },
};
use crate::fields::rational_function_field::RationalFunction;
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Computes `[s]φ` in the additive arithmetic of reduced endomorphisms
    /// by the double-and-add algorithm.
    ///
    /// The returned sum type keeps the additive zero endomorphism explicit and
    /// also preserves the Schoof-relevant branch where one secant/tangent
    /// denominator ceases to be a unit modulo `g(x)`.
    ///
    /// Complexity: `Θ((log s) m^2)`. It performs `Θ(log s)` additions/doublings
    /// for `s = scalar`. Each nontrivial step costs `Θ(m^2)` field operations for
    /// `m = deg g`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn scalar_mul_reduced_endomorphism(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        scalar: u128,
        value: &ReducedEndomorphism<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        if scalar == 0 {
            return ReducedEndomorphismAdditiveResult::Zero;
        }

        let mut remaining = scalar;
        let mut accumulator: Option<ReducedEndomorphism<F>> = None;
        let mut addend = value.clone();

        while remaining > 0 {
            if remaining & 1 == 1 {
                let accumulator_result = accumulator
                    .map(ReducedEndomorphismAdditiveResult::Value)
                    .unwrap_or(ReducedEndomorphismAdditiveResult::Zero);
                let next = accumulator_result.combine_with_value(self, quotient, addend.clone());
                accumulator = match next {
                    Ok(ReducedEndomorphismAdditiveResult::Zero) => None,
                    Ok(ReducedEndomorphismAdditiveResult::Value(value)) => Some(value),
                    Ok(ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd })
                    | Err(witness_gcd) => {
                        return ReducedEndomorphismAdditiveResult::NonUnitDenominator {
                            witness_gcd,
                        };
                    }
                };
            }
            remaining >>= 1;
            if remaining > 0 {
                addend = match self.double_reduced_endomorphism(quotient, &addend) {
                    ReducedEndomorphismAdditiveResult::Value(next) => next,
                    ReducedEndomorphismAdditiveResult::Zero => {
                        remaining = 0;
                        continue;
                    }
                    ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd } => {
                        return ReducedEndomorphismAdditiveResult::NonUnitDenominator {
                            witness_gcd,
                        };
                    }
                };
            }
        }
        accumulator
            .map(ReducedEndomorphismAdditiveResult::Value)
            .unwrap_or(ReducedEndomorphismAdditiveResult::Zero)
    }

    /// Computes `[q] id` for the reduced identity endomorphism, where
    /// `q = #F_q` is the size of the represented finite base field.
    ///
    /// Complexity:  `Θ((log q) m^2)`.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn q_times_reduced_identity_endomorphism(
        &self,
        quotient: &ReducedCurveQuotient<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        let identity = ReducedEndomorphism::identity(quotient);
        let field_order = F::order().expect("finite field order should fit in u128");
        self.scalar_mul_reduced_endomorphism(quotient, field_order, &identity)
    }

    /// Computes `[s] id` on `E[\ell]` by first forming the generic-point
    /// multiple `[s mod \ell](x, y)` in the full short-Weierstrass function
    /// field and only then reducing that affine pair modulo `g(x)`.
    ///
    /// This avoids depending on the reduced affine-addition ladder for scalar
    /// multiples inside the odd-prime Schoof step.
    ///
    /// Complexity: `Θ((log \ell) M)` generic-point operations in `F(E)`, where
    /// `M` is the current function-field arithmetic cost, plus two quotient
    /// reductions modulo `g(x)`.
    pub(crate) fn scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        odd_prime: usize,
        scalar: u128,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        let reduced_scalar = scalar % odd_prime as u128;
        if reduced_scalar == 0 {
            return ReducedEndomorphismAdditiveResult::Zero;
        }

        let function_field = ShortWeierstrassFunctionField::<F>::new(self.clone());
        let generic_multiple = function_field
            .generic_point_multiple(
                u64::try_from(reduced_scalar).expect("odd-prime Schoof scalars should fit in u64"),
            )
            .expect("generic-point scalar multiplication should succeed on valid curves");

        let Some(x_value) = generic_multiple.x() else {
            return ReducedEndomorphismAdditiveResult::Zero;
        };
        let Some(y_value) = generic_multiple.y() else {
            return ReducedEndomorphismAdditiveResult::Zero;
        };

        match self.reduce_generic_affine_pair_to_reduced_endomorphism(quotient, x_value, y_value) {
            Ok(value) => ReducedEndomorphismAdditiveResult::Value(value),
            Err(witness_gcd) => {
                ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd }
            }
        }
    }

    /// Computes `[s]φ` on `E[\ell]` by reducing the canonical generic-point
    /// pullback of `[s mod \ell]` and composing it with `φ`.
    ///
    /// Complexity: the cost of
    /// [`Self::scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion`]
    /// plus one reduced composition `Θ(m^3)` when the scalar multiple is
    /// non-zero.
    pub(crate) fn scalar_multiple_of_reduced_endomorphism_on_odd_torsion(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        odd_prime: usize,
        scalar: u128,
        value: &ReducedEndomorphism<F>,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        match self.scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion(
            quotient, odd_prime, scalar,
        ) {
            ReducedEndomorphismAdditiveResult::Zero => ReducedEndomorphismAdditiveResult::Zero,
            ReducedEndomorphismAdditiveResult::Value(scalar_map) => {
                ReducedEndomorphismAdditiveResult::Value(scalar_map.compose(quotient, value))
            }
            ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd } => {
                ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd }
            }
        }
    }

    fn reduce_generic_affine_pair_to_reduced_endomorphism(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        x_value: &ShortWeierstrassFunction<F>,
        y_value: &ShortWeierstrassFunction<F>,
    ) -> Result<ReducedEndomorphism<F>, DensePolynomial<F>> {
        if !x_value.b_part().is_zero() || !y_value.a_part().is_zero() {
            panic!("generic scalar multiples should stay in the affine shape (a(x), b(x) y)");
        }

        let x_map =
            self.reduce_regular_rational_function_mod_quotient(quotient, x_value.a_part())?;
        let y_scale =
            self.reduce_regular_rational_function_mod_quotient(quotient, y_value.b_part())?;
        Ok(ReducedEndomorphism::new(quotient, x_map, y_scale))
    }

    fn reduce_regular_rational_function_mod_quotient(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        value: &RationalFunction<F>,
    ) -> Result<DensePolynomial<F>, DensePolynomial<F>> {
        match quotient.try_invert_poly(value.denominator()) {
            QuotientInverseResult::Inverse(inverse) => {
                Ok(quotient.reduce_poly(&value.numerator().mul(&inverse)))
            }
            QuotientInverseResult::NonUnit { witness_gcd } => Err(witness_gcd),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        short_weierstrass::division_polynomials::DivisionPolynomialForm,
        short_weierstrass::schoof::{
            ReducedCurveQuotient, ReducedEndomorphism, ReducedEndomorphismAdditiveResult,
        },
    };
    use crate::fields::{Fp, traits::Field};
    use crate::polynomials::DensePolynomial;

    type F7 = Fp<7>;
    type F43 = Fp<43>;

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
    fn zero_scalar_returns_the_additive_zero_endomorphism() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let value = ReducedEndomorphism::identity(&quotient);

        assert_eq!(
            curve.scalar_mul_reduced_endomorphism(&quotient, 0, &value),
            ReducedEndomorphismAdditiveResult::Zero
        );
    }

    #[test]
    fn one_scalar_returns_the_original_endomorphism() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let value = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(6)]),
        );

        assert_eq!(
            curve.scalar_mul_reduced_endomorphism(&quotient, 1, &value),
            ReducedEndomorphismAdditiveResult::Value(value)
        );
    }

    #[test]
    fn two_scalar_matches_the_doubling_helper() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let value = ReducedEndomorphism::identity(&quotient);

        let scalar_mul = curve.scalar_mul_reduced_endomorphism(&quotient, 2, &value);
        let doubled = curve.double_reduced_endomorphism(&quotient, &value);

        match (scalar_mul, doubled) {
            (
                ReducedEndomorphismAdditiveResult::Value(from_scalar_mul),
                ReducedEndomorphismAdditiveResult::Value(from_double),
            ) => {
                assert_eq!(from_scalar_mul, from_double);
            }
            _ => {
                panic!("doubling the chosen identity endomorphism should stay in the affine branch")
            }
        }
    }

    #[test]
    fn two_scalar_of_a_two_torsion_affine_endomorphism_returns_zero() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let value = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::zero()]),
        );

        assert_eq!(
            curve.scalar_mul_reduced_endomorphism(&quotient, 2, &value),
            ReducedEndomorphismAdditiveResult::Zero
        );
    }

    #[test]
    fn q_times_identity_uses_the_same_scalar_multiplication_story() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);

        assert_eq!(
            curve.q_times_reduced_identity_endomorphism(&quotient),
            curve.scalar_mul_reduced_endomorphism(&quotient, 7, &identity)
        );
    }

    #[test]
    fn canonical_two_and_three_times_identity_match_reduced_arithmetic_in_the_seven_torsion_quotient()
     {
        let curve = ShortWeierstrassCurve::<F43>::new(F43::from_i64(-10), F43::from_i64(-10))
            .expect("sample F43 curve should be smooth");
        let DivisionPolynomialForm::InX(psi_seven) = curve
            .division_polynomial(7)
            .expect("psi_7 should exist over F43")
        else {
            panic!("psi_7 should be an x-polynomial for odd ell");
        };
        let quotient = ReducedCurveQuotient::new(curve.clone(), psi_seven)
            .expect("psi_7 should define a reduced quotient");
        let identity = ReducedEndomorphism::identity(&quotient);
        let ReducedEndomorphismAdditiveResult::Value(doubled) =
            curve.scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion(&quotient, 7, 2)
        else {
            panic!("canonical [2]id should stay affine");
        };
        let canonical =
            curve.scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion(&quotient, 7, 3);
        let reduced = curve.add_reduced_endomorphisms(&quotient, &doubled, &identity);
        assert_eq!(canonical, reduced);
    }
}
