use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_ideals::PrimeNormIdeal,
};
use crate::numerics::{
    chinese_remainder::{Congruence, solve_coprime_congruences},
    positive_mod_biguint,
};

/// Direction convention for the prime-ideal to form-class bridge.
///
/// The current bridge stores `PrimeNormIdeal::root_mod_ell()` as a square root
/// `r` of the order discriminant `Δ` modulo `ℓ`. It then chooses the unique
/// `b mod 2ℓ` with `b ≡ r (mod ℓ)` and `b ≡ Δ (mod 2)`, and constructs the
/// primitive form
///
/// (ℓ, b, (b² - Δ)/(4ℓ)).
///
/// Under the standard form-to-ideal convention
/// `(a,b,c) ↦ (a, (-b + √Δ)/2)`, this form records the same local prime above
/// `ℓ` selected by `r`. Conjugating the ideal changes the root to the other
/// root and reduces to the inverse form class.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum IdealFormConvention {
    /// The reduced form represents the supplied prime ideal class.
    RepresentsIdeal,
    /// The reduced form represents the inverse of the supplied ideal class.
    #[allow(dead_code)]
    RepresentsInverseIdeal,
}

/// Form-class label associated to one supported prime-norm ideal.
///
/// This report is intentionally crate-private while the class-group action
/// layer is still staged. It records both the raw prime-norm form with leading
/// coefficient `ℓ` and the reduced representative accepted by the current
/// [`QuadraticClassGroup`] implementation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IdealFormCorrespondence {
    raw_form: BinaryQuadraticForm,
    reduced_form: BinaryQuadraticForm,
    convention: IdealFormConvention,
}

impl IdealFormCorrespondence {
    /// Builds the form-class label attached to `ideal`.
    ///
    /// Let `Δ` be the discriminant of the ideal's order and let `r` be the
    /// stored root of `x² ≡ Δ (mod ℓ)`. The algorithm chooses the unique
    /// representative `b` modulo `2ℓ` satisfying:
    ///
    /// `b ≡ r (mod ℓ)` and `b ≡ Δ (mod 2)`,
    ///
    /// then forms `(ℓ, b, (b² − Δ)/(4ℓ))` and reduces it by Gauss reduction.
    /// The class group is constructed only to validate that the reduced form
    /// lies in the expected discriminant.
    ///
    /// Complexity: `Θ(log ℓ + R + h(D))`, where `R` is the cost of
    /// positive-definite reduction and `h(D)` is the current enumeration-based
    /// validation cost for the class group of discriminant `D`.
    pub(crate) fn from_prime_norm_ideal(
        ideal: &PrimeNormIdeal,
    ) -> Result<Self, BinaryQuadraticFormError> {
        let discriminant = ideal.order().discriminant().value();
        let middle =
            Self::middle_coefficient_for_root(discriminant, ideal.norm(), ideal.root_mod_ell());
        let raw_form = Self::raw_form_for_middle_coefficient(discriminant, ideal.norm(), middle);
        let reduced_form = raw_form.reduce_positive_definite()?;

        let class_group = QuadraticClassGroup::new(ideal.order().discriminant().clone())?;
        if !class_group.contains_reduced_form(&reduced_form) {
            return Err(BinaryQuadraticFormError::NotReducedPositiveDefinite);
        }

        Ok(Self {
            raw_form,
            reduced_form,
            convention: IdealFormConvention::RepresentsIdeal,
        })
    }

    /// Returns the unreduced form `(ℓ,b,c)` produced directly from the root.
    pub(crate) fn raw_form(&self) -> &BinaryQuadraticForm {
        &self.raw_form
    }

    /// Returns the reduced representative of the associated form class.
    pub(crate) fn reduced_form(&self) -> &BinaryQuadraticForm {
        &self.reduced_form
    }

    /// Returns the convention used by this correspondence.
    pub(crate) fn convention(&self) -> IdealFormConvention {
        self.convention
    }

    fn middle_coefficient_for_root(discriminant: &BigInt, ell: &BigUint, root: &BigUint) -> BigInt {
        let discriminant_parity = positive_mod_biguint(discriminant, &BigUint::from(2u8));
        let root_congruence = Congruence::new(root.clone(), ell.clone())
            .expect("prime-norm ideals have positive norm at least 2");
        let parity_congruence = Congruence::new(discriminant_parity, BigUint::from(2u8))
            .expect("mod 2 should define a valid congruence");
        let solution = solve_coprime_congruences(&[root_congruence, parity_congruence])
            .expect("odd prime norm and parity modulus should be coprime");

        BigInt::from(solution.residue().clone())
    }

    fn raw_form_for_middle_coefficient(
        discriminant: &BigInt,
        ell: &BigUint,
        middle: BigInt,
    ) -> BinaryQuadraticForm {
        let leading = BigInt::from(ell.clone());
        BinaryQuadraticForm::from_leading_middle_discriminant(leading, middle, discriminant)
            .expect("selected ideal root should produce an integral prime-norm form")
    }
}
