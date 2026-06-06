use core::{marker::PhantomData, num::NonZeroU32};

use crate::fields::{
    errors::FieldError,
    polynomial_field::PolynomialModulus,
    traits::{Field, FiniteField},
};
use crate::polynomials::{DensePolynomial, PolynomialError};

type BaseElem<S> = <<S as ExtensionFieldSpec>::Base as Field>::Elem;
type DenseTriple<S> = (
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
    DensePolynomial<<S as ExtensionFieldSpec>::Base>,
);

/// Defines an educational quadratic extension over a prime-field family.
///
/// This helper is meant for small examples, tests, and walkthrough code that
/// wants a concrete type-level presentation of
/// `Fp<P>[x] / (x^2 - d)`
/// without rewriting the same [`ExtensionFieldSpec`] boilerplate.
///
/// The generated specification validates the modulus through
/// [`PolynomialModulus::check_field_modulus_requirements`], so the caller is
/// still responsible for choosing a value of `d` that is genuinely
/// non-square when a true quadratic field extension is desired.
///
/// Example:
///
/// ```ignore
/// use elliptic_algorithms_lab::fields::{Field, Fp};
///
/// type F19 = Fp<19>;
///
/// elliptic_algorithms_lab::fields::define_fp_quadratic_extension!(
///     spec: F19Sqrt2Spec,
///     field: F19Sqrt2,
///     base: F19,
///     non_residue: 2,
///     name: "F19(sqrt(2))",
/// );
/// ```
#[macro_export]
macro_rules! define_fp_quadratic_extension {
    (
        spec: $spec:ident,
        field: $field:ident,
        base: $base:ty,
        non_residue: $non_residue:expr,
        name: $name:expr $(,)?
    ) => {
        struct $spec;

        impl $crate::fields::ExtensionFieldSpec for $spec {
            type Base = $base;

            const NAME: &'static str = $name;

            fn defining_modulus() -> $crate::fields::PolynomialModulus<Self::Base> {
                $crate::fields::PolynomialModulus::<Self::Base>::new(vec![
                    <Self::Base as $crate::fields::Field>::from_i64(-($non_residue)),
                    <Self::Base as $crate::fields::Field>::zero(),
                    <Self::Base as $crate::fields::Field>::one(),
                ])
                .expect("x^2 - d should be a valid structural modulus")
            }

            fn check_field_conditions() -> Result<(), $crate::fields::FieldError> {
                Self::defining_modulus().check_field_modulus_requirements()
            }
        }

        type $field = $crate::fields::ExtensionField<$spec>;
    };
}

/// Defines an educational quadratic extension over `Q`.
///
/// This helper is meant for small exact examples such as `Q(sqrt(2))` that
/// would otherwise repeat the same [`ExtensionFieldSpec`] boilerplate.
///
/// The generated specification validates the modulus through
/// [`PolynomialModulus::check_field_modulus_requirements`], so it remains
/// honest about whether `x^2 - d` really defines a field extension.
///
/// Example:
///
/// ```ignore
/// elliptic_algorithms_lab::fields::define_q_quadratic_extension!(
///     spec: QSqrt2Spec,
///     field: QSqrt2,
///     radicand: 2,
///     name: "Q(sqrt(2))",
/// );
/// ```
#[macro_export]
macro_rules! define_q_quadratic_extension {
    (
        spec: $spec:ident,
        field: $field:ident,
        radicand: $radicand:expr,
        name: $name:expr $(,)?
    ) => {
        struct $spec;

        impl $crate::fields::ExtensionFieldSpec for $spec {
            type Base = $crate::fields::Q;

            const NAME: &'static str = $name;

            fn defining_modulus() -> $crate::fields::PolynomialModulus<Self::Base> {
                $crate::fields::PolynomialModulus::<Self::Base>::new(vec![
                    <$crate::fields::Q as $crate::fields::Field>::from_i64(-($radicand)),
                    <$crate::fields::Q as $crate::fields::Field>::zero(),
                    <$crate::fields::Q as $crate::fields::Field>::one(),
                ])
                .expect("x^2 - d should be a valid structural modulus")
            }

            fn check_field_conditions() -> Result<(), $crate::fields::FieldError> {
                Self::defining_modulus().check_field_modulus_requirements()
            }
        }

        type $field = $crate::fields::ExtensionField<$spec>;
    };
}

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

    fn dense_modulus() -> DensePolynomial<S::Base> {
        DensePolynomial::<S::Base>::new(Self::modulus().coefficients().to_vec())
    }

    fn element_to_dense(element: &ExtensionFieldElement<S>) -> DensePolynomial<S::Base> {
        DensePolynomial::<S::Base>::new(element.coefficients.clone())
    }

    fn from_dense_reduced(polynomial: DensePolynomial<S::Base>) -> ExtensionFieldElement<S> {
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

/// Canonical representative of an element in `Base[x] / (m(x))`.
///
/// The stored coefficients are kept in ascending degree order and are trimmed
/// to remove trailing zero coefficients. They do not carry a runtime field
/// descriptor because the ambient quotient is already determined by the type
/// parameter `S`.
pub struct ExtensionFieldElement<S: ExtensionFieldSpec> {
    coefficients: Vec<BaseElem<S>>,
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

impl<S: ExtensionFieldSpec> Field for ExtensionField<S> {
    const IS_ALGEBRAICALLY_CLOSED: bool = S::IS_ALGEBRAICALLY_CLOSED;

    type Elem = ExtensionFieldElement<S>;

    fn characteristic() -> u64 {
        S::Base::characteristic()
    }

    fn zero() -> Self::Elem {
        Self::zero_element()
    }

    fn one() -> Self::Elem {
        Self::one_element()
    }

    fn from_i64(n: i64) -> Self::Elem {
        Self::from_base(S::Base::from_i64(n))
    }

    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::add_elements(x, y)
    }

    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::sub_elements(x, y)
    }

    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Self::mul_elements(x, y)
    }

    fn neg(x: &Self::Elem) -> Self::Elem {
        Self::neg_element(x)
    }

    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        Self::inverse_element(x).ok()
    }

    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        let left = Self::reduce(x);
        let right = Self::reduce(y);

        left.coefficients.len() == right.coefficients.len()
            && left
                .coefficients
                .iter()
                .zip(&right.coefficients)
                .all(|(lhs, rhs)| S::Base::eq(lhs, rhs))
    }

    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::inverse_element(x)
    }

    fn elem_from_u64(value: u64) -> Self::Elem {
        Self::from_base(S::Base::elem_from_u64(value))
    }
}

impl<S> FiniteField for ExtensionField<S>
where
    S: ExtensionFieldSpec,
    S::Base: FiniteField,
{
    fn extension_degree() -> NonZeroU32 {
        let base_degree = <S::Base as FiniteField>::extension_degree().get();
        let step_degree = Self::extension_degree().get();
        let total_degree = base_degree
            .checked_mul(step_degree)
            .expect("finite-field extension degree should fit in u32");
        NonZeroU32::new(total_degree).expect("finite-field extension degree must stay non-zero")
    }

    fn check_structure() -> Result<(), FieldError> {
        <S::Base as FiniteField>::check_structure()?;
        Self::check_structure()
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use proptest::prelude::*;

    use crate::fields::{ComplexApprox, Field, FieldError, FiniteField, Fp, PolynomialModulus, Q};
    use crate::fields::{ExtensionField, ExtensionFieldElement, ExtensionFieldSpec};
    use crate::proptest_support::{
        ProptestF17Sqrt3Field, ProptestF17TowerField, tower_element_case,
    };

    type F17 = Fp<17>;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    crate::fields::define_q_quadratic_extension!(
        spec: QSqrt2Spec,
        field: QSqrt2,
        radicand: 2,
        name: "Q(sqrt(2))",
    );

    struct QSqrt2ISpec;

    impl ExtensionFieldSpec for QSqrt2ISpec {
        type Base = QSqrt2;

        const NAME: &'static str = "Q(sqrt(2), i)";

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<QSqrt2>::new(vec![QSqrt2::one(), QSqrt2::zero(), QSqrt2::one()])
                .expect("x^2 + 1 should be a valid structural modulus")
        }

        fn check_field_conditions() -> Result<(), FieldError> {
            // The crate does not yet expose a generic irreducibility backend
            // over algebraic extension bases, so this tower step is accepted
            // from the known mathematics of Q(sqrt(2), i).
            Ok(())
        }
    }

    type QSqrt2I = ExtensionField<QSqrt2ISpec>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F17Sqrt3Spec,
        field: F17Sqrt3,
        base: F17,
        non_residue: 3,
        name: "F17(sqrt(3))",
    );

    struct F17TowerSpec;

    impl ExtensionFieldSpec for F17TowerSpec {
        type Base = F17Sqrt3;

        const NAME: &'static str = "F17(sqrt(3))(u)";

        fn defining_modulus() -> PolynomialModulus<Self::Base> {
            PolynomialModulus::<F17Sqrt3>::new(vec![
                F17Sqrt3::one(),
                F17Sqrt3::one(),
                F17Sqrt3::zero(),
                F17Sqrt3::one(),
            ])
            .expect("tower modulus should be structurally valid")
        }

        fn check_field_conditions() -> Result<(), FieldError> {
            Ok(())
        }
    }

    type F17Tower = ExtensionField<F17TowerSpec>;

    #[test]
    fn extension_field_spec_can_be_constructed_as_zero_sized_field() {
        let _field = QSqrt2::new().expect("q(sqrt(2)) should validate");
        assert_eq!(QSqrt2::name(), "Q(sqrt(2))");
        assert_eq!(QSqrt2::extension_degree().get(), 2);
    }

    #[test]
    fn extension_field_reduces_x_squared_to_base_relation() {
        let x_squared = QSqrt2::element(vec![Q::zero(), Q::zero(), Q::one()]);
        let reduced = QSqrt2::reduce(&x_squared);

        assert_eq!(reduced.coefficients().len(), 1);
        assert!(Q::eq(&reduced.coefficients()[0], &q(2, 1)));
    }

    #[test]
    fn extension_field_addition_and_multiplication_work_over_q_sqrt2() {
        let one_plus_x = QSqrt2::element(vec![Q::one(), Q::one()]);
        let three_minus_x = QSqrt2::element(vec![Q::from_i64(3), Q::from_i64(-1)]);

        let sum = QSqrt2::add(&one_plus_x, &three_minus_x);
        let product = QSqrt2::mul(&one_plus_x, &three_minus_x);

        assert_eq!(sum.coefficients().len(), 1);
        assert!(Q::eq(&sum.coefficients()[0], &q(4, 1)));

        assert_eq!(product.coefficients().len(), 2);
        assert!(Q::eq(&product.coefficients()[0], &q(1, 1)));
        assert!(Q::eq(&product.coefficients()[1], &q(2, 1)));
    }

    #[test]
    fn extension_field_inverse_is_computed_via_polynomial_euclid() {
        let element = QSqrt2::element(vec![Q::one(), Q::one()]);
        let inverse = QSqrt2::inverse(&element).expect("1 + sqrt(2) should be invertible");
        let check = QSqrt2::mul(&element, &inverse);

        assert!(QSqrt2::eq(&check, &QSqrt2::one()));
        assert_eq!(inverse.coefficients().len(), 2);
        assert!(Q::eq(&inverse.coefficients()[0], &q(-1, 1)));
        assert!(Q::eq(&inverse.coefficients()[1], &q(1, 1)));
    }

    #[test]
    fn extension_field_elements_compare_by_quotient_equivalence() {
        let x_squared =
            ExtensionFieldElement::<QSqrt2Spec>::new(vec![Q::zero(), Q::zero(), Q::one()]);
        let two = QSqrt2::element(vec![Q::from_i64(2)]);

        assert_eq!(x_squared, two);
    }

    #[test]
    fn extension_field_can_be_used_as_base_of_another_extension() {
        let i = QSqrt2I::element(vec![QSqrt2::zero(), QSqrt2::one()]);
        let i_squared = QSqrt2I::mul(&i, &i);
        let minus_one = QSqrt2I::element(vec![QSqrt2::from_i64(-1)]);

        assert_eq!(i_squared, minus_one);
    }

    #[test]
    fn extension_field_embeds_base_elements_in_towers() {
        let embedded = QSqrt2I::from_i64(2);
        let expected = QSqrt2I::element(vec![QSqrt2::from_i64(2)]);

        assert_eq!(embedded, expected);
    }

    #[test]
    fn finite_field_metadata_multiplies_through_extension_towers() {
        assert_eq!(F17Sqrt3::characteristic(), 17);
        assert_eq!(F17Sqrt3::extension_degree().get(), 2);
        assert_eq!(<F17Tower as FiniteField>::extension_degree().get(), 6);
    }

    #[test]
    fn extension_field_structure_rejects_constant_modulus_after_normalization() {
        struct BadSpec;

        impl ExtensionFieldSpec for BadSpec {
            type Base = F17;

            fn defining_modulus() -> PolynomialModulus<Self::Base> {
                PolynomialModulus::<F17>::new(vec![F17::one(), F17::zero()])
                    .expect("structural constructor should accept the raw coefficients")
            }
        }

        let error = ExtensionField::<BadSpec>::check_structure()
            .expect_err("effectively constant modulus should be rejected");
        assert!(matches!(error, FieldError::InvalidPolynomialModulus));
    }

    #[test]
    fn extension_field_non_invertible_zero_is_rejected() {
        let error = QSqrt2::inverse(&QSqrt2::zero()).expect_err("zero should not be invertible");
        assert!(matches!(error, FieldError::DivisionByZero));
    }

    #[test]
    fn extension_field_algebraic_closedness_is_forwarded_from_the_spec() {
        fn algebraically_closed<F: Field>() -> bool {
            F::IS_ALGEBRAICALLY_CLOSED
        }

        struct ComplexLinearSpec;

        impl ExtensionFieldSpec for ComplexLinearSpec {
            type Base = ComplexApprox;

            const NAME: &'static str = "C(a)";
            const IS_ALGEBRAICALLY_CLOSED: bool = true;

            fn defining_modulus() -> PolynomialModulus<Self::Base> {
                PolynomialModulus::<ComplexApprox>::new(vec![
                    ComplexApprox::from_i64(-1),
                    ComplexApprox::one(),
                ])
                .expect("linear modulus should be structurally valid")
            }
        }

        assert!(algebraically_closed::<ExtensionField<ComplexLinearSpec>>());
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(28))]

        #[test]
        fn property_extension_field_reduction_is_idempotent(
            coefficients in prop::collection::vec(0u64..17, 0..=5),
        ) {
            let element = F17Sqrt3::element(
                coefficients
                    .into_iter()
                    .map(F17::elem_from_u64)
                    .collect(),
            );

            let reduced_once = F17Sqrt3::reduce(&element);
            let reduced_twice = F17Sqrt3::reduce(&reduced_once);

            prop_assert_eq!(reduced_once, reduced_twice);
        }

        #[test]
        fn property_extension_field_nonzero_elements_are_invertible(
            coefficients in prop::collection::vec(0u64..17, 0..=2),
        ) {
            let element = F17Sqrt3::element(
                coefficients
                    .into_iter()
                    .map(F17::elem_from_u64)
                    .collect(),
            );
            let reduced = F17Sqrt3::reduce(&element);

            prop_assume!(!F17Sqrt3::is_zero(&reduced));

            let inverse = F17Sqrt3::inverse(&reduced).expect("non-zero extension-field element should be invertible");
            let one = F17Sqrt3::mul(&reduced, &inverse);

            prop_assert!(F17Sqrt3::eq(&one, &F17Sqrt3::one()));
        }

        #[test]
        fn property_base_embedding_into_quadratic_extension_is_a_field_homomorphism(
            case in tower_element_case(),
        ) {
            let embedded_left = ProptestF17Sqrt3Field::from_base(case.base_left);
            let embedded_right = ProptestF17Sqrt3Field::from_base(case.base_right);
            let embedded_sum = ProptestF17Sqrt3Field::add(&embedded_left, &embedded_right);
            let embedded_product = ProptestF17Sqrt3Field::mul(&embedded_left, &embedded_right);

            prop_assert_eq!(
                embedded_sum,
                ProptestF17Sqrt3Field::from_base(F17::add(&case.base_left, &case.base_right))
            );
            prop_assert_eq!(
                embedded_product,
                ProptestF17Sqrt3Field::from_base(F17::mul(&case.base_left, &case.base_right))
            );
        }

        #[test]
        fn property_quadratic_base_embedding_into_tower_is_a_field_homomorphism(
            case in tower_element_case(),
        ) {
            let embedded_left = ProptestF17TowerField::from_base(case.quadratic_left.clone());
            let embedded_right = ProptestF17TowerField::from_base(case.quadratic_right.clone());
            let embedded_sum = ProptestF17TowerField::add(&embedded_left, &embedded_right);
            let embedded_product = ProptestF17TowerField::mul(&embedded_left, &embedded_right);

            prop_assert_eq!(
                embedded_sum,
                ProptestF17TowerField::from_base(ProptestF17Sqrt3Field::add(
                    &case.quadratic_left,
                    &case.quadratic_right,
                ))
            );
            prop_assert_eq!(
                embedded_product,
                ProptestF17TowerField::from_base(ProptestF17Sqrt3Field::mul(
                    &case.quadratic_left,
                    &case.quadratic_right,
                ))
            );
        }

        #[test]
        fn property_semantic_tower_elements_stay_canonically_reduced_under_arithmetic(
            case in tower_element_case(),
        ) {
            let sum = ProptestF17TowerField::add(&case.tower_left, &case.tower_right);
            let product = ProptestF17TowerField::mul(&case.tower_left, &case.tower_right);

            prop_assert_eq!(ProptestF17TowerField::reduce(&sum), sum);
            prop_assert_eq!(ProptestF17TowerField::reduce(&product), product);
        }
    }
}
