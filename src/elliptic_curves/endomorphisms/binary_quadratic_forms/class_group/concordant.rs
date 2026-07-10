use num_bigint::BigInt;
use num_traits::Zero;

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};
use crate::numerics::{
    chinese_remainder::{ChineseRemainderSolution, Congruence, combine_compatible_congruences},
    gcd_bigint, positive_bigint_to_biguint, positive_mod_biguint,
};

impl QuadraticClassGroup {
    /// Composes two concordant representatives.
    ///
    /// The input forms `(a, b, c)` and `(a', b', c')` are concordant when
    /// `gcd(a, a', (b + b')/2) = 1`.
    ///
    /// In this case Dirichlet composition first forms `A = aa'`, chooses the
    /// middle coefficient `B` from the congruences `B ≡ b (mod 2a)`,
    /// `B ≡ b′ (mod 2a')`, and `B² ≡ D (mod 4A)`. Then `C = (B² − D)/(4A)`,
    /// and the resulting form `(A, B, C)` is returned after Gauss reduction.
    ///
    /// The inputs must be primitive positive-definite representatives of this
    /// class group, but they need not already be reduced. The returned form is
    /// reduced in the positive-definite Gauss convention.
    ///
    /// Complexity: one compatible CRT combination, followed by `Θ(k)` candidate
    /// middle coefficients, where `k` is the number of residue lifts inspected
    /// before `B² − D` is divisible by `4aa′`.
    pub(crate) fn compose_concordant_forms(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
    ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError> {
        self.validate_member(left)?;
        self.validate_member(right)?;
        debug_assert!(self.are_concordant(left, right));

        let leading = left.a() * right.a();
        let middle_residue = self.concordant_middle_residue(left, right)?;
        let middle = self
            .find_integral_middle_coefficient(&leading, middle_residue)
            .expect(
                "valid concordant primitive inputs should admit an integral middle coefficient",
            );
        BinaryQuadraticForm::from_leading_middle_discriminant(
            leading,
            middle,
            self.discriminant().value(),
        )
        .expect("chosen concordant middle coefficient should make the constant integral")
        .reduce_positive_definite()
    }

    fn are_concordant(&self, left: &BinaryQuadraticForm, right: &BinaryQuadraticForm) -> bool {
        let shared_leading_gcd = gcd_bigint(left.a(), right.a());
        let middle_sum_half = (left.b() + right.b()) / BigInt::from(2u8);
        let concordance_gcd = gcd_bigint(&shared_leading_gcd, &middle_sum_half);

        concordance_gcd == BigInt::from(1u8)
    }

    fn concordant_middle_residue(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
    ) -> Result<ChineseRemainderSolution, BinaryQuadraticFormError> {
        let left_modulus = positive_bigint_to_biguint(&(BigInt::from(2u8) * left.a()));
        let right_modulus = positive_bigint_to_biguint(&(BigInt::from(2u8) * right.a()));
        let left_residue = positive_mod_biguint(left.b(), &left_modulus);
        let right_residue = positive_mod_biguint(right.b(), &right_modulus);
        let left_solution = ChineseRemainderSolution::new(left_residue, left_modulus);
        let right_congruence = Congruence::new(right_residue, right_modulus)
            .expect("validated positive-definite forms have positive leading coefficients");

        Ok(combine_compatible_congruences(&left_solution, &right_congruence).expect(
            "validated concordant forms with the same discriminant have compatible middle congruences",
        ))
    }

    fn find_integral_middle_coefficient(
        &self,
        leading: &BigInt,
        residue: ChineseRemainderSolution,
    ) -> Option<BigInt> {
        let period = positive_bigint_to_biguint(&(BigInt::from(2u8) * leading));
        let mut candidate = residue.residue().clone();
        let denominator = BigInt::from(4u8) * leading;

        while candidate < period {
            let middle = BigInt::from(candidate.clone());
            let numerator = &middle * &middle - self.discriminant().value();
            if (&numerator % &denominator).is_zero() {
                return Some(middle);
            }
            candidate += residue.modulus();
        }
        None
    }
}
