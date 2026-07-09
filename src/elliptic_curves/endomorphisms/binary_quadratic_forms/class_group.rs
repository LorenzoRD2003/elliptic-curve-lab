use num_bigint::BigInt;
use num_traits::{One, Signed};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError},
    quadratic_orders::QuadraticDiscriminant,
};

/// Enumerative scaffold for the class group of an imaginary quadratic order.
///
/// This first layer stores only a negative quadratic-order discriminant `D`
/// and enumerates the primitive reduced positive-definite forms of
/// discriminant `D`. It does not yet implement composition of classes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuadraticClassGroup {
    discriminant: QuadraticDiscriminant,
}

impl QuadraticClassGroup {
    /// Builds the class-group scaffold for one imaginary quadratic discriminant.
    ///
    /// The first enumerator supports exactly the discriminants of imaginary
    /// quadratic orders: `D < 0` and `D ≡ 0, 1 (mod 4)`.
    pub fn new(discriminant: QuadraticDiscriminant) -> Result<Self, BinaryQuadraticFormError> {
        if !discriminant.is_negative() {
            return Err(BinaryQuadraticFormError::NotNegativeDiscriminant);
        }

        if !discriminant.is_congruent_to_0_mod_4() && !discriminant.is_congruent_to_1_mod_4() {
            return Err(BinaryQuadraticFormError::NotQuadraticOrderDiscriminant);
        }

        Ok(Self { discriminant })
    }

    /// Returns the negative quadratic-order discriminant.
    pub fn discriminant(&self) -> &QuadraticDiscriminant {
        &self.discriminant
    }

    /// Enumerates primitive reduced positive-definite forms of discriminant `D`.
    ///
    /// For a reduced positive-definite form `(a,b,c)` of discriminant `D < 0`,
    /// the inequalities imply `1 ≤ a` and `3a² ≤ |D|`. The enumeration uses
    /// that finite bound, tests all `−a ≤ b ≤ a`, and keeps the primitive
    /// forms satisfying `b² - 4ac = D` together with the reduced convention
    /// implemented by [`BinaryQuadraticForm::is_reduced_positive_definite`].
    pub fn enumerate_reduced_forms(&self) -> Vec<BinaryQuadraticForm> {
        let discriminant = self.discriminant.value();
        let abs_discriminant = discriminant.abs();
        let mut forms = Vec::new();
        let mut a = BigInt::one();

        while BigInt::from(3u8) * &a * &a <= abs_discriminant {
            let mut b = -&a;
            while b <= a {
                if let Some(form) = self.reduced_form_with_coefficients(&a, &b) {
                    forms.push(form);
                }
                b += 1u8;
            }
            a += 1u8;
        }
        forms
    }

    fn reduced_form_with_coefficients(
        &self,
        a: &BigInt,
        b: &BigInt,
    ) -> Option<BinaryQuadraticForm> {
        let discriminant = self.discriminant.value();
        let numerator = b * b - discriminant;
        let denominator = BigInt::from(4u8) * a;

        if (&numerator % &denominator) != BigInt::from(0u8) {
            return None;
        }

        let c = numerator / denominator;
        let form = BinaryQuadraticForm::new(a.clone(), b.clone(), c);

        (form.discriminant() == *discriminant
            && form.is_primitive()
            && form.is_reduced_positive_definite())
        .then_some(form)
    }
}
