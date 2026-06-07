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
