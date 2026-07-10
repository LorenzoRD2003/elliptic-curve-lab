use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};
use crate::numerics::{
    chinese_remainder::{ChineseRemainderSolution, Congruence, combine_compatible_congruences},
    gcd_bigint, positive_bigint_to_biguint, positive_mod_biguint,
};

impl QuadraticClassGroup {
    /// Composes two reduced representatives in this class group.
    ///
    /// The inputs must be primitive reduced positive-definite binary quadratic
    /// forms of this group's discriminant `D`. The result is the reduced
    /// representative of the product class in the proper class group of the
    /// imaginary quadratic order of discriminant `D`.
    ///
    /// For reduced primitive positive-definite forms `(a, b, c)` and
    /// `(a', b', c')` of the same discriminant `D`, set
    ///
    /// `g = gcd(a, a′, (b + b′)/2)` and `A = aa′ / g²`.
    ///
    /// The middle coefficient `B` is chosen so that
    /// `B ≡ b  (mod 2a/g)`, `B ≡ b′ (mod 2a′/g)`, and `B² ≡ D (mod 4A)`.
    ///
    /// Then `C = (B² − D)/(4A)`, and the resulting form `(A, B, C)` is
    /// returned after positive-definite Gauss reduction. This is the classical
    /// Dirichlet/Gauss composition route.
    ///
    /// Roadmap: this method intentionally exposes only the group operation,
    /// not the chosen composition engine. A later milestone may replace or
    /// supplement the current classical Dirichlet/Gauss engine with a NUCOMP
    /// strategy while preserving this austere caller-facing API.
    ///
    /// Complexity: exact gcd computations on the leading coefficients and
    /// middle sum, one compatible CRT combination, and `Θ(k)` inspected lifts
    /// of the CRT residue, where `k` is bounded by the natural period `2A`,
    /// plus positive-definite Gauss reduction of the product representative.
    pub fn compose(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
    ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError> {
        self.validate_reduced_member(left)?;
        self.validate_reduced_member(right)?;

        let shared_factor = self.dirichlet_common_factor(left, right);
        if shared_factor.is_one() {
            return self.compose_concordant_forms(left, right);
        }

        let leading = (left.a() * right.a()) / (&shared_factor * &shared_factor);
        let middle_residue = self.dirichlet_middle_residue(left, right, &shared_factor);
        let middle = self
            .find_dirichlet_middle_coefficient(&leading, middle_residue)
            .expect("valid reduced primitive inputs should admit a Dirichlet middle coefficient");
        BinaryQuadraticForm::from_leading_middle_discriminant(
            leading,
            middle,
            self.discriminant().value(),
        )
        .expect("chosen Dirichlet middle coefficient should make the constant integral")
        .reduce_positive_definite()
    }

    fn dirichlet_common_factor(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
    ) -> BigInt {
        let shared_leading_gcd = gcd_bigint(left.a(), right.a());
        let middle_sum_half = (left.b() + right.b()) / BigInt::from(2u8);

        gcd_bigint(&shared_leading_gcd, &middle_sum_half)
    }

    fn dirichlet_middle_residue(
        &self,
        left: &BinaryQuadraticForm,
        right: &BinaryQuadraticForm,
        common_factor: &BigInt,
    ) -> ChineseRemainderSolution {
        let left_modulus =
            positive_bigint_to_biguint(&((BigInt::from(2u8) * left.a()) / common_factor));
        let right_modulus =
            positive_bigint_to_biguint(&((BigInt::from(2u8) * right.a()) / common_factor));
        let left_residue = positive_mod_biguint(left.b(), &left_modulus);
        let right_residue = positive_mod_biguint(right.b(), &right_modulus);
        let left_solution = ChineseRemainderSolution::new(left_residue, left_modulus);
        let right_congruence = Congruence::new(right_residue, right_modulus)
            .expect("validated reduced positive-definite forms have positive leading coefficients");

        combine_compatible_congruences(&left_solution, &right_congruence).expect(
            "validated reduced forms with the same discriminant have compatible Dirichlet middle congruences",
        )
    }

    fn find_dirichlet_middle_coefficient(
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
