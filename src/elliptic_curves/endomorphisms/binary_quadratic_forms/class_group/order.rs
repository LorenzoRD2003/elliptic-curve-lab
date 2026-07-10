use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};

impl QuadraticClassGroup {
    /// Returns the order of a reduced form class in this quadratic class group.
    ///
    /// The input must already be the reduced positive-definite representative
    /// used by this group. The algorithm starts from the principal class and
    /// repeatedly composes by `form` until it returns to the principal class.
    /// Thus the result is the least `n ≥ 1` such that `[form]^n = 1`.
    ///
    /// Complexity: `O(h(D) · C)`, where `h(D)` is the number of reduced forms
    /// of discriminant `D` and `C` is the cost of one Dirichlet composition
    /// followed by Gauss reduction. This is intentionally simple and suited to
    /// the small educational discriminants currently used by the repository.
    pub fn order_of_reduced_form(
        &self,
        form: &BinaryQuadraticForm,
    ) -> Result<usize, BinaryQuadraticFormError> {
        self.generated_subgroup(form)
            .map(|subgroup| subgroup.order())
    }

    pub(super) fn principal_reduced_form(
        &self,
    ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError> {
        let middle = if self.discriminant().is_congruent_to_0_mod_4() {
            BigInt::zero()
        } else {
            BigInt::one()
        };

        BinaryQuadraticForm::from_leading_middle_discriminant(
            BigInt::one(),
            middle,
            self.discriminant().value(),
        )
        .ok_or(BinaryQuadraticFormError::ClassOrderNotFound)?
        .reduce_positive_definite()
    }
}
