use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::BinaryQuadraticFormError, quadratic_orders::QuadraticDiscriminant,
};

/// Enumerative scaffold for the class group of an imaginary quadratic order.
///
/// This first layer stores only a negative quadratic-order discriminant `D`
/// and enumerates the primitive reduced positive-definite forms of
/// discriminant `D`. It does not yet expose composition of classes.
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
}
