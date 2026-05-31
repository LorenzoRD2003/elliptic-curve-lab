use core::marker::PhantomData;

use crate::fields::{errors::FieldError, traits::Field};

/// Modulus polynomial used to define a quotient of `F[x]`.
#[derive(Clone, Debug)]
pub struct PolynomialModulus<F: Field> {
    coefficients: Vec<F::Elem>,
    field: PhantomData<F>,
}

impl<F: Field> PolynomialModulus<F> {
    /// Creates a modulus polynomial from coefficients in ascending order.
    ///
    /// This constructor currently performs only the weakest structural check:
    /// the modulus must be non-constant.
    ///
    /// That is enough to talk about a quotient algebra `F[x] / (m(x))`, but it
    /// is not enough to guarantee that the quotient behaves as a field.
    ///
    /// For quotient-field constructions, the modulus polynomial is typically
    /// expected to be irreducible over the coefficient field `F`. That stronger
    /// requirement is intentionally not enforced here yet, because:
    ///
    /// - irreducibility testing has not been implemented
    /// - this type is also useful for generic quotient-algebra scaffolding
    /// - the project still wants to keep the low-level constructor lightweight
    ///
    /// In other words:
    ///
    /// - `new(...)` means “this is a syntactically valid non-constant modulus”
    /// - it does **not** yet mean “this polynomial is suitable for a true
    ///   quotient field”
    pub fn new(coefficients: Vec<F::Elem>) -> Result<Self, FieldError> {
        if coefficients.len() < 2 {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self {
            coefficients,
            field: PhantomData,
        })
    }

    /// Returns the raw coefficients.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }

    /// Returns the degree of the modulus polynomial.
    pub fn degree(&self) -> usize {
        self.coefficients.len().saturating_sub(1)
    }

    /// Checks whether the modulus satisfies the stronger conditions usually
    /// required for quotient-field constructions.
    ///
    /// Once implemented, this method should be the natural place to ask
    /// questions such as irreducibility over the coefficient field `F`.
    pub fn check_field_modulus_requirements(&self) -> Result<(), FieldError> {
        todo!("irreducibility checking for polynomial moduli is intentionally deferred")
    }
}

/// Element of a quotient algebra `F[x] / (m(x))`.
///
/// Mathematically, this type is meant to model a polynomial representative
/// together with the modulus polynomial that defines the quotient relation.
/// Two representatives should be understood as denoting the same quotient
/// element when they differ by a multiple of `m(x)`.
///
/// In other words, if `m(x)` is the chosen modulus polynomial, then
/// polynomials `f(x)` and `g(x)` represent the same quotient element whenever:
///
/// `f(x) - g(x) = q(x) * m(x)`
///
/// for some polynomial `q(x)` over the same coefficient field `F`.
///
/// This struct deliberately stores:
///
/// - a raw representative polynomial through `coefficients`
/// - the modulus polynomial that defines the quotient space
///
/// The current representation uses coefficients in ascending degree order.
/// For example:
///
/// - `[a0, a1, a2]` represents `a0 + a1*x + a2*x^2`
/// - a modulus `[b0, b1, b2]` represents `b0 + b1*x + b2*x^2`
///
/// This type is still scaffold-oriented:
///
/// - representatives are not automatically reduced yet
/// - equivalence between unreduced representatives is not normalized yet
/// - no irreducibility checks are enforced yet
///
/// That means `PolynomialFieldElement<F>` should currently be read as
/// “a future quotient-field element with an explicit representative” rather
/// than as a fully normalized algebraic value.
#[derive(Clone, Debug)]
pub struct PolynomialFieldElement<F: Field> {
    coefficients: Vec<F::Elem>,
    modulus: PolynomialModulus<F>,
}

impl<F: Field> PolynomialFieldElement<F> {
    /// Builds a quotient element from a polynomial representative.
    ///
    /// The `coefficients` slice stores the chosen representative in ascending
    /// degree order, while `modulus` defines the quotient relation.
    ///
    /// At the moment this constructor validates only the basic structural
    /// requirement that the modulus polynomial is non-constant. It does not yet:
    ///
    /// - reduce the representative modulo the modulus
    /// - remove trailing zero coefficients
    /// - verify that the modulus is irreducible
    ///
    /// Those behaviors are intentionally deferred until the polynomial
    /// arithmetic layer is implemented in a more complete way.
    pub fn new(
        coefficients: Vec<F::Elem>,
        modulus: PolynomialModulus<F>,
    ) -> Result<Self, FieldError> {
        if modulus.coefficients().is_empty() {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        Ok(Self {
            coefficients,
            modulus,
        })
    }

    /// Returns the representative coefficients.
    ///
    /// The coefficients are stored in ascending degree order. The returned
    /// slice exposes the current explicit representative, which may still be
    /// unreduced in this scaffold phase.
    pub fn coefficients(&self) -> &[F::Elem] {
        &self.coefficients
    }

    /// Returns the defining modulus.
    ///
    /// Conceptually, this is the polynomial `m(x)` in the quotient
    /// construction `F[x] / (m(x))`.
    pub fn modulus(&self) -> &PolynomialModulus<F> {
        &self.modulus
    }

    /// Reduces the representative polynomial modulo the field relation.
    ///
    /// Once implemented, this method should turn the current representative
    /// into a canonical remainder modulo the defining polynomial.
    pub fn reduce(&self) -> Result<Self, FieldError> {
        todo!("quotient reduction will be implemented together with polynomial division")
    }

    /// Checks the modulus for the extra conditions required by a field.
    ///
    /// In a true quotient field, the modulus polynomial must satisfy stronger
    /// conditions than in a generic quotient algebra, most notably
    /// irreducibility over the coefficient field.
    pub fn check_field_conditions(&self) -> Result<(), FieldError> {
        self.modulus.check_field_modulus_requirements()
    }
}

#[cfg(test)]
mod tests {
    use super::{PolynomialFieldElement, PolynomialModulus};
    use crate::fields::{Field, Fp};

    type F17 = Fp<17>;

    #[test]
    fn quotient_modulus_requires_degree_at_least_one() {
        let error = PolynomialModulus::<F17>::new(vec![F17::elem_from_u64(42)])
            .expect_err("constant modulus should fail");
        assert!(matches!(
            error,
            crate::fields::FieldError::InvalidPolynomialModulus
        ));
    }

    #[test]
    fn quotient_modulus_preserves_coefficients_and_degree() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(3),
            F17::elem_from_u64(0),
            F17::elem_from_u64(5),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        let coefficients = modulus.coefficients();
        assert_eq!(coefficients.len(), 4);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(0)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(5)));
        assert!(F17::eq(&coefficients[3], &F17::elem_from_u64(1)));
        assert_eq!(modulus.degree(), 3);
    }

    #[test]
    fn quotient_element_preserves_representative_coefficients() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(9),
                F17::elem_from_u64(4),
                F17::elem_from_u64(15),
            ],
            modulus,
        )
        .expect("element should exist");

        let coefficients = element.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(9)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(4)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(15)));
    }

    #[test]
    fn quotient_element_exposes_defining_modulus() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(2),
            F17::elem_from_u64(3),
            F17::elem_from_u64(4),
        ])
        .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(vec![F17::elem_from_u64(8)], modulus)
            .expect("element should exist");

        let stored_modulus = element.modulus();
        assert_eq!(stored_modulus.degree(), 2);
        let coefficients = stored_modulus.coefficients();
        assert_eq!(coefficients.len(), 3);
        assert!(F17::eq(&coefficients[0], &F17::elem_from_u64(2)));
        assert!(F17::eq(&coefficients[1], &F17::elem_from_u64(3)));
        assert!(F17::eq(&coefficients[2], &F17::elem_from_u64(4)));
    }

    #[test]
    fn quotient_element_can_store_zero_representative() {
        let modulus =
            PolynomialModulus::<F17>::new(vec![F17::elem_from_u64(1), F17::elem_from_u64(1)])
                .expect("modulus should exist");

        let element = PolynomialFieldElement::<F17>::new(Vec::new(), modulus)
            .expect("zero representative should be allowed");

        assert!(element.coefficients().is_empty());
    }

    #[test]
    #[ignore = "field-modulus checks are intentionally scaffold-only"]
    fn quotient_modulus_field_requirements_placeholder() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        modulus
            .check_field_modulus_requirements()
            .expect("placeholder");
    }

    #[test]
    #[ignore = "irreducibility checking is intentionally scaffold-only"]
    fn quotient_element_field_conditions_placeholder() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(vec![F17::elem_from_u64(5)], modulus)
            .expect("element should exist");
        element.check_field_conditions().expect("placeholder");
    }

    #[test]
    #[ignore = "quotient-field arithmetic is scaffold-only"]
    fn quotient_reduction_placeholder() {
        let modulus = PolynomialModulus::<F17>::new(vec![
            F17::elem_from_u64(1),
            F17::elem_from_u64(0),
            F17::elem_from_u64(1),
        ])
        .expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(
            vec![
                F17::elem_from_u64(1),
                F17::elem_from_u64(2),
                F17::elem_from_u64(3),
            ],
            modulus,
        )
        .expect("element should exist");
        let _ = element.reduce().expect("placeholder");
    }
}
