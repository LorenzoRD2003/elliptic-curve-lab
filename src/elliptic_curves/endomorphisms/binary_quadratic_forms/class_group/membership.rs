use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};

impl QuadraticClassGroup {
    /// Validates a reduced representative of this imaginary quadratic class group.
    ///
    /// A representative belongs to the class-group when its discriminant is the
    /// group discriminant `D`, it is primitive, and it is already in the reduced
    /// positive-definite convention used by Gauss reduction.
    pub(super) fn validate_reduced_member(
        &self,
        form: &BinaryQuadraticForm,
    ) -> Result<(), BinaryQuadraticFormError> {
        self.validate_member(form)?;

        if !form.is_reduced_positive_definite() {
            return Err(BinaryQuadraticFormError::NotReducedPositiveDefinite);
        }
        Ok(())
    }

    /// Validates a primitive positive-definite representative of this class group.
    ///
    /// This accepts non-reduced representatives, which appear naturally after
    /// applying a proper equivalence before Dirichlet composition.
    pub(super) fn validate_member(
        &self,
        form: &BinaryQuadraticForm,
    ) -> Result<(), BinaryQuadraticFormError> {
        if form.discriminant() != *self.discriminant().value() {
            return Err(BinaryQuadraticFormError::ClassGroupDiscriminantMismatch);
        } else if !form.is_primitive() {
            return Err(BinaryQuadraticFormError::NotPrimitive);
        } else if !form.is_positive_definite() {
            return Err(BinaryQuadraticFormError::NotPositiveDefinite);
        }
        Ok(())
    }
}
