use num_bigint::BigInt;
use num_traits::{One, Signed};

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, QuadraticClassGroup,
};

impl QuadraticClassGroup {
    /// Enumerates primitive reduced positive-definite forms of discriminant `D`.
    ///
    /// For a reduced positive-definite form `(a,b,c)` of discriminant `D < 0`,
    /// the inequalities imply `1 ≤ a` and `3a² ≤ |D|`. The enumeration uses
    /// that finite bound, tests all `−a ≤ b ≤ a`, and keeps the primitive
    /// forms satisfying `b² - 4ac = D` together with the reduced convention
    /// implemented by [`BinaryQuadraticForm::is_reduced_positive_definite`].
    pub fn enumerate_reduced_forms(&self) -> Vec<BinaryQuadraticForm> {
        let discriminant = self.discriminant().value();
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
        let discriminant = self.discriminant().value();
        let numerator = b * b - discriminant;
        let denominator = BigInt::from(4u8) * a;

        if (&numerator % &denominator) != BigInt::from(0u8) {
            return None;
        }

        let c = numerator / denominator;
        let form = BinaryQuadraticForm::new(a.clone(), b.clone(), c);

        self.validate_reduced_member(&form).ok().map(|()| form)
    }
}
