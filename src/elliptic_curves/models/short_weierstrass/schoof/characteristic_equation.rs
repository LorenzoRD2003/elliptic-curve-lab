use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::schoof::{
        ReducedCurveQuotient, ReducedEndomorphism, ReducedEndomorphismAdditiveResult,
    },
};
use crate::fields::traits::FiniteField;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Forms the reduced Schoof characteristic-equation candidate
    /// `π^2 - [c]π + [q] id`.
    ///
    /// Here `π` is represented by the supplied reduced endomorphism
    /// `frobenius`, `c` is one candidate trace residue, and `[q] id` denotes
    /// the `q`-fold additive multiple of the identity endomorphism in the
    /// reduced additive arithmetic.
    ///
    /// Complexity: `Θ((log c + log q + m) m^2)` field operations, where
    /// `m = deg g`, dominated by scalar multiplications in the additive
    /// arithmetic and one reduced composition `π^2`.
    pub(crate) fn reduced_characteristic_equation_candidate(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        odd_prime: usize,
        frobenius: &ReducedEndomorphism<F>,
        frobenius_squared: &ReducedEndomorphism<F>,
        q_term: &ReducedEndomorphismAdditiveResult<F>,
        candidate_trace: u128,
    ) -> ReducedEndomorphismAdditiveResult<F> {
        let trace_term = self.scalar_multiple_of_reduced_endomorphism_on_odd_torsion(
            quotient,
            odd_prime,
            candidate_trace,
            frobenius,
        );
        let negated_trace_term = trace_term.additive_inverse();

        let after_subtraction = match ReducedEndomorphismAdditiveResult::Value(
            frobenius_squared.clone(),
        )
        .combine(self, quotient, negated_trace_term)
        {
            Ok(result) => result,
            Err(witness_gcd) => {
                return ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd };
            }
        };

        match after_subtraction.combine(self, quotient, q_term.clone()) {
            Ok(result) => result,
            Err(witness_gcd) => {
                ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd }
            }
        }
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
            DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient")
    }

    #[test]
    fn characteristic_equation_candidate_can_return_the_additive_zero() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);
        let frobenius_squared = identity.compose(&quotient, &identity);
        let q_term = curve.q_times_reduced_identity_endomorphism(&quotient);

        assert_eq!(
            curve.reduced_characteristic_equation_candidate(
                &quotient,
                17,
                &identity,
                &frobenius_squared,
                &q_term,
                8,
            ),
            ReducedEndomorphismAdditiveResult::Zero
        );
    }

    #[test]
    fn characteristic_equation_candidate_can_return_a_nonzero_value() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);
        let frobenius_squared = identity.compose(&quotient, &identity);
        let q_term = curve.q_times_reduced_identity_endomorphism(&quotient);

        assert_eq!(
            curve.reduced_characteristic_equation_candidate(
                &quotient,
                17,
                &identity,
                &frobenius_squared,
                &q_term,
                7,
            ),
            ReducedEndomorphismAdditiveResult::Value(identity)
        );
    }

    #[test]
    fn characteristic_equation_candidate_can_return_a_higher_scalar_multiple() {
        let curve = sample_curve();
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);
        let frobenius_squared = identity.compose(&quotient, &identity);
        let q_term = curve.q_times_reduced_identity_endomorphism(&quotient);
        let expected = curve.scalar_mul_reduced_endomorphism(&quotient, 8, &identity);

        assert_eq!(
            curve.reduced_characteristic_equation_candidate(
                &quotient,
                17,
                &identity,
                &frobenius_squared,
                &q_term,
                0,
            ),
            expected
        );
    }

    #[test]
    fn characteristic_equation_candidate_propagates_non_unit_witnesses() {
        let curve = sample_curve();
        let quotient = ReducedCurveQuotient::new(
            sample_curve(),
            DensePolynomial::new(vec![F7::from_i64(-1), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient");
        let identity = ReducedEndomorphism::identity(&quotient);
        let frobenius_squared = identity.compose(&quotient, &identity);
        let q_term = curve.q_times_reduced_identity_endomorphism(&quotient);

        let ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd } =
            curve.reduced_characteristic_equation_candidate(
                &quotient,
                17,
                &identity,
                &frobenius_squared,
                &q_term,
                2,
            )
        else {
            panic!("the chosen quotient should hit the non-unit tangent branch");
        };

        assert_eq!(
            witness_gcd,
            DensePolynomial::new(vec![F7::one(), F7::one()])
        );
    }
}
