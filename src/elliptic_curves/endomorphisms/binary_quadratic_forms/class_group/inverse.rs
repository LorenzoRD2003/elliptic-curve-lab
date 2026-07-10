use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
};

impl QuadraticClassGroup {
    /// Returns the reduced representative of the inverse class.
    ///
    /// For a primitive positive-definite form `(a, b, c)`, the inverse proper
    /// class is represented by the conjugate form `(a, −b, c)`. This method
    /// validates that the input is a reduced representative of this class group
    /// and then reduces that conjugate representative.
    ///
    /// Complexity: one reduced-membership validation plus positive-definite
    /// Gauss reduction of the conjugate form.
    pub fn inverse(
        &self,
        form: &BinaryQuadraticForm,
    ) -> Result<BinaryQuadraticForm, BinaryQuadraticFormError> {
        self.validate_reduced_member(form)?;
        form.conjugate().reduce_positive_definite()
    }
}
