use core::num::NonZeroU32;

use crate::fields::{errors::FieldError, polynomial_field::PolynomialModulus, traits::Field};
use crate::polynomials::IrreducibilityBackend;

/// Metadata for a field extension presented as a quotient `F[x] / (m(x))`.
///
/// This descriptor is intentionally generic over the base field `F`. An
/// extension field does not need to be an extension of a finite field:
/// infinite fields such as `Q` also admit non-trivial algebraic extensions.
///
/// The quotient should be understood as a true field only when `m(x)` satisfies
/// the stronger conditions required for that construction, most notably
/// irreducibility over `F`. The current crate now exposes those checks through
/// [`ExtensionFieldDescriptor::check_field_conditions`] when the base field has
/// an irreducibility backend.
#[derive(Clone, Debug)]
pub struct ExtensionFieldDescriptor<F: Field> {
    /// Defining polynomial for the quotient presentation `F[x] / (m(x))`.
    pub modulus: PolynomialModulus<F>,
}

impl<F: Field> ExtensionFieldDescriptor<F> {
    /// Builds an extension descriptor from its defining modulus polynomial.
    ///
    /// The constructor currently performs only the structural check already
    /// enforced by [`PolynomialModulus::new`]: the modulus must be
    /// non-constant.
    ///
    /// This is enough to describe a quotient of the polynomial ring `F[x]`.
    /// It is not, by itself, enough to certify that the quotient is a genuine
    /// field. Use [`ExtensionFieldDescriptor::check_field_conditions`] for the
    /// stronger validation step when the base field supports it.
    pub fn new(modulus: PolynomialModulus<F>) -> Result<Self, FieldError> {
        if modulus.coefficients().len() < 2 {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self { modulus })
    }

    /// Returns the defining modulus polynomial.
    pub fn modulus(&self) -> &PolynomialModulus<F> {
        &self.modulus
    }

    /// Returns the algebraic degree of the extension over the base field `F`.
    ///
    /// For a quotient presented as `F[x] / (m(x))`, this degree is the degree
    /// of the defining polynomial `m(x)`.
    pub fn extension_degree(&self) -> NonZeroU32 {
        let degree = self.modulus.degree() as u32;
        NonZeroU32::new(degree)
            .expect("extension-field descriptors require a non-constant modulus polynomial")
    }

    /// Checks whether the defining modulus satisfies the extra conditions
    /// normally required for a true field extension.
    ///
    /// This method is available only when the base field implements the
    /// polynomial irreducibility backend used by the crate.
    pub fn check_field_conditions(&self) -> Result<(), FieldError>
    where
        F: IrreducibilityBackend,
    {
        self.modulus.check_field_modulus_requirements()
    }
}

/// Operational view of an extension field presented as `F[x] / (m(x))`.
///
/// This type is the natural home for operations that depend on the defining
/// quotient relation. It plays the same conceptual role that `Fp<P>` plays for
/// prime fields:
///
/// - `ExtensionField<F>` stores the ambient field description
/// - `ExtensionFieldElement<F>` stores only a representative value
///
/// This separation keeps element values lightweight while avoiding APIs that
/// require callers to thread a raw descriptor through every operation.
#[derive(Clone, Debug)]
pub struct ExtensionField<F: Field> {
    /// Descriptor of the ambient extension field.
    pub descriptor: ExtensionFieldDescriptor<F>,
}

impl<F: Field> ExtensionField<F> {
    /// Builds an operational extension field from a descriptor.
    pub fn new(descriptor: ExtensionFieldDescriptor<F>) -> Result<Self, FieldError> {
        if descriptor.modulus().coefficients().len() < 2 {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self { descriptor })
    }

    /// Builds an operational extension field directly from a modulus
    /// polynomial.
    pub fn from_modulus(modulus: PolynomialModulus<F>) -> Result<Self, FieldError> {
        Self::new(ExtensionFieldDescriptor::new(modulus)?)
    }

    /// Returns the descriptor of the ambient extension.
    pub fn descriptor(&self) -> &ExtensionFieldDescriptor<F> {
        &self.descriptor
    }

    /// Returns the algebraic degree of the extension over the base field.
    pub fn extension_degree(&self) -> NonZeroU32 {
        self.descriptor.extension_degree()
    }

    /// Builds an element of this extension from a polynomial representative.
    ///
    /// The representative is stored exactly as provided. Reduction modulo the
    /// defining polynomial is deferred to the arithmetic layer.
    pub fn element(
        &self,
        coefficients: Vec<F::Elem>,
    ) -> Result<ExtensionFieldElement<F>, FieldError> {
        ExtensionFieldElement::new(coefficients)
    }

    /// Adds two elements in this extension.
    pub fn add(
        &self,
        _left: &ExtensionFieldElement<F>,
        _right: &ExtensionFieldElement<F>,
    ) -> Result<ExtensionFieldElement<F>, FieldError> {
        todo!(
            "coordinate-wise addition will be implemented once the polynomial-field arithmetic is in place"
        )
    }

    /// Multiplies two elements in this extension before quotient reduction.
    pub fn mul(
        &self,
        _left: &ExtensionFieldElement<F>,
        _right: &ExtensionFieldElement<F>,
    ) -> Result<ExtensionFieldElement<F>, FieldError> {
        todo!(
            "extension-field multiplication is deferred until polynomial multiplication APIs stabilize"
        )
    }

    /// Reduces an element representative modulo the defining polynomial.
    pub fn reduce(
        &self,
        _element: &ExtensionFieldElement<F>,
    ) -> Result<ExtensionFieldElement<F>, FieldError> {
        todo!("polynomial reduction strategy is intentionally left open in the scaffold")
    }

    /// Computes the multiplicative inverse when the element is invertible.
    pub fn inverse(
        &self,
        _element: &ExtensionFieldElement<F>,
    ) -> Result<ExtensionFieldElement<F>, FieldError> {
        todo!("extension-field inversion will depend on the chosen Euclidean polynomial backend")
    }
}

/// Element of a quotient presentation `F[x] / (m(x))`.
///
/// The base field is any `F: Field`, not necessarily finite. This makes the
/// type appropriate both for classic finite-field extensions and for algebraic
/// extensions of infinite fields such as the rationals.
///
/// The struct stores only the polynomial representative in ascending degree
/// order.
///
/// For example:
///
/// - `coefficients = [a0, a1]` means `a0 + a1*x`
/// - `coefficients = [a0, a1, a2]` means `a0 + a1*x + a2*x^2`
///
/// The ambient extension metadata is intentionally not stored inside the
/// element. This mirrors the separation already used by prime-field elements
/// in the crate.
#[derive(Clone, Debug)]
pub struct ExtensionFieldElement<F: Field> {
    /// Coefficients of the stored polynomial representative in ascending degree
    /// order.
    pub coefficients: Vec<F::Elem>,
}

impl<F: Field> ExtensionFieldElement<F> {
    /// Builds an element from a polynomial representative.
    pub fn new(coefficients: Vec<F::Elem>) -> Result<Self, FieldError> {
        Ok(Self { coefficients })
    }

    /// Returns the stored polynomial representative.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::{ExtensionField, ExtensionFieldDescriptor, ExtensionFieldElement};
    use crate::fields::{Field, Fp, PolynomialModulus, Q};

    type F17 = Fp<17>;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn extension_descriptor_preserves_modulus_and_degree() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        let descriptor =
            ExtensionFieldDescriptor::<F17>::new(modulus).expect("descriptor should exist");

        assert_eq!(descriptor.extension_degree().get(), 2);
        let coefficients = descriptor.modulus().coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(1)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(0)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(1)));
    }

    #[test]
    fn extension_field_exposes_descriptor_and_degree() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let field = ExtensionField::<F17>::from_modulus(modulus).expect("field should exist");

        assert_eq!(field.extension_degree().get(), 2);
        let stored_modulus = field.descriptor().modulus().coefficients();
        assert_eq!(stored_modulus.len(), 3);
        assert!(F17::eq(&stored_modulus[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&stored_modulus[1], &F17::elem_from_u64(0)));
        assert!(F17::eq(&stored_modulus[2], &F17::elem_from_u64(1)));
    }

    #[test]
    fn extension_element_preserves_representative_coefficients() {
        let element =
            ExtensionFieldElement::<F17>::new(vec![F17::elem_from_u64(8), F17::elem_from_u64(5)])
                .expect("element should exist");

        let coefficients = element.coefficients();
        assert_eq!(coefficients.len(), 2);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(8)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(5)));
    }

    #[test]
    fn extension_element_can_store_zero_representative() {
        let element =
            ExtensionFieldElement::<F17>::new(Vec::new()).expect("zero representative is valid");

        assert!(element.coefficients().is_empty());
    }

    #[test]
    fn extension_field_builds_elements_without_storing_descriptor_in_them() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let field = ExtensionField::<F17>::from_modulus(modulus).expect("field should exist");

        let element = field
            .element(vec![F17::elem_from_u64(9), F17::elem_from_u64(4)])
            .expect("element should exist");

        let coefficients = element.coefficients();
        assert_eq!(coefficients.len(), 2);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(9)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(4)));
    }

    #[test]
    fn extension_descriptor_is_generic_over_infinite_base_fields_too() {
        let modulus = PolynomialModulus::<Q>::new(vec![Q::from_i64(-2), Q::zero(), Q::one()])
            .expect("modulus should exist");

        let descriptor = ExtensionFieldDescriptor::<Q>::new(modulus)
            .expect("descriptor should exist over the rationals too");

        assert_eq!(descriptor.extension_degree().get(), 2);
        let coefficients = descriptor.modulus().coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(Q::eq(&coefficients[0], &q(-2, 1)));
        assert!(Q::eq(&coefficients[1], &q(0, 1)));
        assert!(Q::eq(&coefficients[2], &q(1, 1)));
    }

    #[test]
    fn extension_field_is_generic_over_infinite_base_fields_too() {
        let modulus = PolynomialModulus::<Q>::new(vec![Q::from_i64(-2), Q::zero(), Q::one()])
            .expect("modulus should exist");
        let field = ExtensionField::<Q>::from_modulus(modulus).expect("field should exist");

        assert_eq!(field.extension_degree().get(), 2);
        let stored_modulus = field.descriptor().modulus().coefficients();
        assert_eq!(stored_modulus.len(), 3);
        assert!(Q::eq(&stored_modulus[0], &q(-2, 1)));
        assert!(Q::eq(&stored_modulus[1], &q(0, 1)));
        assert!(Q::eq(&stored_modulus[2], &q(1, 1)));
    }

    #[test]
    fn extension_field_can_validate_true_field_conditions() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let descriptor =
            ExtensionFieldDescriptor::<F17>::new(modulus).expect("descriptor should exist");
        let field = ExtensionField::<F17>::new(descriptor).expect("field should exist");

        field
            .descriptor()
            .check_field_conditions()
            .expect("irreducible modulus should pass");
    }

    #[test]
    #[ignore = "extension arithmetic is intentionally deferred in this scaffold"]
    fn extension_field_addition_is_available_once_polynomial_arithmetic_exists() {
        let modulus =
            PolynomialModulus::<F17>::new(vec![F17::elem_from_u64(1), F17::elem_from_u64(1)])
                .expect("modulus should exist");
        let field = ExtensionField::<F17>::from_modulus(modulus).expect("field should exist");

        let left =
            ExtensionFieldElement::<F17>::new(vec![F17::elem_from_u64(1)]).expect("left element");
        let right = left.clone();

        let _sum = field
            .add(&left, &right)
            .expect("addition should eventually work");
    }
}
