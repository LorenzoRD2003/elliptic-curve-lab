use crate::fields::traits::*;
use crate::fields::{error::FieldError, polynomial_field::PolynomialModulus};
use crate::polynomials::DensePolynomial;

/// Static specification of an algebraic field extension presented as
/// `Base[x] / (m(x))`.
///
/// This trait is the key to making extension fields fit into the same
/// type-level `Field` API as prime fields, rationals, and complex numbers.
///
/// Instead of storing the defining modulus as runtime state inside each field
/// object or element value, the extension data is attached to a dedicated
/// specification type. That makes the quotient presentation part of the type
/// itself, which in turn enables:
///
/// - `impl Field for ExtensionField<S>`
/// - polynomial arithmetic over extension fields without introducing a second
///   runtime-only field trait
/// - towers such as `Q(sqrt(2))` and `Q(sqrt(2), i)` by using one extension
///   field as the base of another
///
/// The contract for implementors is:
///
/// - [`ExtensionFieldSpec::defining_modulus`] must return the polynomial that
///   defines the quotient
/// - [`ExtensionFieldSpec::check_field_conditions`] should validate that the
///   quotient really is a field whenever the current crate has enough
///   infrastructure to decide that question
/// - when the crate does not yet have a generic backend for the base field,
///   this method is also the place to record a mathematically justified manual
///   acceptance of the extension
///
/// In particular, towers over non-prime bases are now possible even before the
/// crate has full irreducibility machinery over arbitrary algebraic extension
/// fields.
pub trait ExtensionFieldSpec {
    /// Base field over which the defining polynomial is written.
    type Base: Field;

    /// Short educational label for the extension family.
    ///
    /// This is used only by documentation and visualization helpers.
    const NAME: &'static str = "unnamed extension field";

    /// Whether the modeled extension field family is algebraically closed.
    ///
    /// Most algebraic extensions in this crate will leave this as `false`.
    /// The hook exists so the specification can state the mathematical truth
    /// explicitly when a degenerate or unusual quotient presentation models an
    /// algebraically closed field family.
    const IS_ALGEBRAICALLY_CLOSED: bool = false;

    /// Returns the defining modulus polynomial `m(x)`.
    fn defining_modulus() -> PolynomialModulus<Self::Base>;

    /// Checks the stronger conditions required for the quotient to behave as a
    /// true field.
    ///
    /// The default implementation performs only the weakest structural check:
    /// the defining polynomial must remain non-constant after dense
    /// normalization.
    ///
    /// When the base field implements a polynomial irreducibility backend, a
    /// typical implementation will delegate to
    /// [`PolynomialModulus::check_field_modulus_requirements`].
    ///
    /// For tower steps whose base field does not yet have a generic
    /// irreducibility backend, this hook can currently be used to document and
    /// accept a mathematically known valid extension manually.
    fn check_field_conditions() -> Result<(), FieldError> {
        let modulus = Self::defining_modulus();
        let dense = DensePolynomial::<Self::Base>::new(modulus.coefficients().to_vec());

        if dense.degree().is_some_and(|degree| degree >= 1) {
            Ok(())
        } else {
            Err(FieldError::InvalidPolynomialModulus)
        }
    }
}
