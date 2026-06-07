use core::{marker::PhantomData, num::NonZeroU32};

use crate::fields::{
    Field,
    errors::FieldError,
    extension_field::{BaseElem, DenseTriple},
    polynomial_field::PolynomialModulus,
};
use crate::polynomials::{DensePolynomial, PolynomialError};

use super::{ExtensionFieldElement, ExtensionFieldSpec};

/// Type-level field backend for an algebraic extension presented as
/// `Base[x] / (m(x))`.
///
/// The field itself carries no runtime data. All ambient information lives in
/// the [`ExtensionFieldSpec`] implementation `S`.
///
/// This makes `ExtensionField<S>` conceptually parallel to `Fp<P>`:
///
/// - `Fp<P>` is a compile-time field family with a static modulus `P`
/// - `ExtensionField<S>` is a compile-time field family with a static defining
///   polynomial supplied by `S`
///
/// Because the extension family is known at the type level, it can itself
/// serve as the base field of another extension. That is the mechanism that
/// enables towers such as:
///
/// - `Q(sqrt(2))`
/// - `Q(sqrt(2), i)`
/// - finite-field towers like `Fp -> Fp2 -> Fp6 -> Fp12`
#[derive(Clone, Copy, Debug)]
pub struct ExtensionField<S: ExtensionFieldSpec>(PhantomData<S>);

impl<S: ExtensionFieldSpec> ExtensionField<S> {
    /// Returns the educational name of the extension family.
    pub fn name() -> &'static str {
        S::NAME
    }

    /// Builds a zero-sized field handle after validating the defining
    /// polynomial.
    ///
    /// This constructor is mostly useful for APIs that prefer to pass around a
    /// field handle explicitly, such as visualization helpers. Arithmetic does
    /// not require an instance because [`Field`] methods remain type-level.
    pub fn new() -> Result<Self, FieldError> {
        Self::check_structure()?;
        Ok(Self(PhantomData))
    }

    /// Returns the defining modulus polynomial of the quotient presentation.
    pub fn modulus() -> PolynomialModulus<S::Base> {
        S::defining_modulus()
    }

    /// Returns the algebraic degree of the extension over its immediate base
    /// field.
    pub fn extension_degree() -> NonZeroU32 {
        let dense_modulus = Self::dense_modulus();
        let degree = dense_modulus
            .degree()
            .expect("extension-field modulus must be non-constant");
        let degree = u32::try_from(degree).expect("extension-field degree should fit in u32");
        NonZeroU32::new(degree).expect("extension-field degree must be non-zero")
    }

    /// Validates the quotient-field conditions recorded by the specification.
    pub fn check_structure() -> Result<(), FieldError> {
        let dense_modulus = Self::dense_modulus();
        if dense_modulus.degree().is_none_or(|degree| degree == 0) {
            return Err(FieldError::InvalidPolynomialModulus);
        }

        S::check_field_conditions()
    }

    /// Builds a quotient representative from coefficients in ascending degree
    /// order and reduces it modulo the defining polynomial.
    pub fn element(coefficients: Vec<BaseElem<S>>) -> ExtensionFieldElement<S> {
        Self::from_dense_reduced(DensePolynomial::<S::Base>::new(coefficients))
    }

    /// Embeds an element of the base field into the extension as a constant
    /// polynomial class.
    pub fn from_base(value: BaseElem<S>) -> ExtensionFieldElement<S> {
        ExtensionFieldElement::new(vec![value])
    }

    /// Returns the additive identity as a canonical quotient representative.
    pub fn zero_element() -> ExtensionFieldElement<S> {
        ExtensionFieldElement::new(Vec::new())
    }

    /// Returns the multiplicative identity as a canonical quotient
    /// representative.
    pub fn one_element() -> ExtensionFieldElement<S> {
        ExtensionFieldElement::new(vec![S::Base::one()])
    }

    /// Returns a canonical representative of the quotient class.
    pub fn reduce(element: &ExtensionFieldElement<S>) -> ExtensionFieldElement<S> {
        Self::from_dense_reduced(Self::element_to_dense(element))
    }

    /// Adds two extension-field elements and reduces the result canonically.
    pub fn add_elements(
        left: &ExtensionFieldElement<S>,
        right: &ExtensionFieldElement<S>,
    ) -> ExtensionFieldElement<S> {
        let sum = Self::element_to_dense(left).add(&Self::element_to_dense(right));
        Self::from_dense_reduced(sum)
    }

    /// Negates an extension-field element coefficient-wise.
    pub fn neg_element(element: &ExtensionFieldElement<S>) -> ExtensionFieldElement<S> {
        ExtensionFieldElement::new(element.coefficients.iter().map(S::Base::neg).collect())
    }

    /// Subtracts two extension-field elements and reduces the result
    /// canonically.
    pub fn sub_elements(
        left: &ExtensionFieldElement<S>,
        right: &ExtensionFieldElement<S>,
    ) -> ExtensionFieldElement<S> {
        Self::add_elements(left, &Self::neg_element(right))
    }

    /// Multiplies two extension-field elements and reduces the product modulo
    /// the defining polynomial.
    pub fn mul_elements(
        left: &ExtensionFieldElement<S>,
        right: &ExtensionFieldElement<S>,
    ) -> ExtensionFieldElement<S> {
        let product = Self::element_to_dense(left).mul(&Self::element_to_dense(right));
        Self::from_dense_reduced(product)
    }

    /// Computes the multiplicative inverse of an element when it exists.
    ///
    /// The current implementation uses the extended Euclidean algorithm in the
    /// polynomial ring `Base[x]`.
    pub fn inverse_element(
        element: &ExtensionFieldElement<S>,
    ) -> Result<ExtensionFieldElement<S>, FieldError> {
        let reduced = Self::reduce(element);
        if reduced.is_zero() {
            return Err(FieldError::DivisionByZero);
        }

        let modulus = Self::dense_modulus();
        let representative = Self::element_to_dense(&reduced);
        let (gcd, bezout, _) =
            Self::extended_gcd(representative, modulus).map_err(Self::map_polynomial_error)?;

        let Some(unit) = gcd.constant_term() else {
            return Err(FieldError::NonInvertibleElement);
        };

        if gcd.degree().is_some_and(|degree| degree > 0) {
            return Err(FieldError::NonInvertibleElement);
        }

        let unit_inverse = S::Base::inverse(unit)?;
        let inverse = bezout.scale(&unit_inverse);
        Ok(Self::from_dense_reduced(inverse))
    }

    pub(super) fn dense_modulus() -> DensePolynomial<S::Base> {
        DensePolynomial::<S::Base>::new(Self::modulus().coefficients().to_vec())
    }

    pub(super) fn element_to_dense(element: &ExtensionFieldElement<S>) -> DensePolynomial<S::Base> {
        DensePolynomial::<S::Base>::new(element.coefficients.clone())
    }

    pub(super) fn from_dense_reduced(
        polynomial: DensePolynomial<S::Base>,
    ) -> ExtensionFieldElement<S> {
        let remainder = polynomial
            .rem(&Self::dense_modulus())
            .expect("extension-field modulus must be a non-zero polynomial");
        ExtensionFieldElement::new(remainder.coefficients().to_vec())
    }

    fn extended_gcd(
        left: DensePolynomial<S::Base>,
        right: DensePolynomial<S::Base>,
    ) -> Result<DenseTriple<S>, PolynomialError> {
        let mut old_r = left;
        let mut r = right;
        let mut old_s = DensePolynomial::<S::Base>::constant(S::Base::one());
        let mut s = DensePolynomial::<S::Base>::new(Vec::new());
        let mut old_t = DensePolynomial::<S::Base>::new(Vec::new());
        let mut t = DensePolynomial::<S::Base>::constant(S::Base::one());

        while !r.is_zero() {
            let (quotient, remainder) = old_r.div_rem(&r)?;
            old_r = r;
            r = remainder;

            let next_s = old_s.sub(&quotient.mul(&s));
            old_s = s;
            s = next_s;

            let next_t = old_t.sub(&quotient.mul(&t));
            old_t = t;
            t = next_t;
        }

        Ok((old_r, old_s, old_t))
    }

    fn map_polynomial_error(error: PolynomialError) -> FieldError {
        match error {
            PolynomialError::DivisionByZeroPolynomial => FieldError::InvalidPolynomialModulus,
            PolynomialError::NonInvertibleLeadingCoefficient => FieldError::NonInvertibleElement,
            _ => FieldError::Unsupported(
                "unexpected polynomial-domain error during extension-field arithmetic",
            ),
        }
    }
}
