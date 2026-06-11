use core::hash::{Hash, Hasher};

use crate::fields::{Field, extension_field::BaseElem};

use super::{ExtensionField, ExtensionFieldSpec};

/// Canonical representative of an element in `Base[x] / (m(x))`.
///
/// The stored coefficients are kept in ascending degree order and are trimmed
/// to remove trailing zero coefficients. They do not carry a runtime field
/// descriptor because the ambient quotient is already determined by the type
/// parameter `S`.
pub struct ExtensionFieldElement<S: ExtensionFieldSpec> {
    pub(super) coefficients: Vec<BaseElem<S>>,
}

impl<S: ExtensionFieldSpec> Clone for ExtensionFieldElement<S> {
    fn clone(&self) -> Self {
        Self {
            coefficients: self.coefficients.clone(),
        }
    }
}

impl<S: ExtensionFieldSpec> core::fmt::Debug for ExtensionFieldElement<S> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ExtensionFieldElement")
            .field("coefficients", &self.coefficients)
            .finish()
    }
}

impl<S: ExtensionFieldSpec> PartialEq for ExtensionFieldElement<S> {
    fn eq(&self, other: &Self) -> bool {
        <ExtensionField<S> as Field>::eq(self, other)
    }
}

impl<S> Eq for ExtensionFieldElement<S>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: Eq,
{
}

impl<S> Hash for ExtensionFieldElement<S>
where
    S: ExtensionFieldSpec,
    BaseElem<S>: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.coefficients.hash(state);
    }
}

impl<S: ExtensionFieldSpec> ExtensionFieldElement<S> {
    /// Builds an element representative from ascending-degree coefficients.
    ///
    /// The constructor trims trailing zeros but does not perform quotient
    /// reduction on its own. Use [`ExtensionField::element`] when a canonical
    /// reduced representative is desired immediately.
    pub fn new(coefficients: Vec<BaseElem<S>>) -> Self {
        Self {
            coefficients: Self::trim_trailing_zeros(coefficients),
        }
    }

    /// Returns the stored coefficients in ascending degree order.
    pub fn coefficients(&self) -> &[BaseElem<S>] {
        &self.coefficients
    }

    /// Returns the degree of the stored representative, if it is non-zero.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }

    /// Returns whether the stored representative is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }

    fn trim_trailing_zeros(mut coefficients: Vec<BaseElem<S>>) -> Vec<BaseElem<S>> {
        while coefficients.last().is_some_and(S::Base::is_zero) {
            coefficients.pop();
        }

        coefficients
    }
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use super::ExtensionFieldElement;
    use crate::fields::Field;
    use crate::fields::extension_field::ExtensionFieldSpec;
    use crate::fields::{Fp, PolynomialModulus};

    type F5 = Fp<5>;

    struct F5Sqrt2Spec;

    impl ExtensionFieldSpec for F5Sqrt2Spec {
        type Base = F5;

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<Self::Base>::new(vec![
                Self::Base::from_i64(-2),
                Self::Base::zero(),
                Self::Base::one(),
            ])
            .expect("x^2 - 2 should be a valid structural modulus over F5")
        }
    }

    fn hash_of<T: Hash>(value: &T) -> u64 {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[test]
    fn trimming_trailing_zeros_preserves_hash() {
        let trimmed = ExtensionFieldElement::<F5Sqrt2Spec>::new(vec![F5::one()]);
        let with_trailing_zero =
            ExtensionFieldElement::<F5Sqrt2Spec>::new(vec![F5::one(), F5::zero()]);

        assert_eq!(trimmed, with_trailing_zero);
        assert_eq!(hash_of(&trimmed), hash_of(&with_trailing_zero));
    }
}
