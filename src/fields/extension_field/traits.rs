use core::num::NonZeroU32;

use crate::fields::{
    FieldError,
    extension_field::{BaseElem, ExtensionField, ExtensionFieldElement, ExtensionFieldSpec},
    traits::{
        CbrtField, EnumerableFiniteField, Field, FiniteField, QuadraticCharacterFiniteField,
        SqrtField,
    },
};

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

impl<S: ExtensionFieldSpec> FiniteField for ExtensionField<S>
where
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

impl<S: ExtensionFieldSpec> EnumerableFiniteField for ExtensionField<S>
where
    S::Base: EnumerableFiniteField,
{
    /// Returns every extension-field element through canonical coefficient
    /// tuples of degree strictly less than the defining modulus degree.
    ///
    /// If the extension is presented as `Base[x] / (m(x))` with
    /// `deg(m) = d`, then every element admits a unique representative
    ///
    /// `a_0 + a_1 x + ... + a_{d-1} x^{d-1}`
    ///
    /// with `a_i in Base`. This implementation enumerates those tuples
    /// directly in deterministic coefficient-lexicographic order.
    fn elements() -> Vec<Self::Elem> {
        let coefficient_slots = ExtensionField::<S>::extension_degree().get() as usize;
        let base_elements = S::Base::elements();
        let total = <Self as FiniteField>::cardinality()
            .and_then(|value| usize::try_from(value).ok())
            .expect("enumerable extension field cardinality should fit in usize");

        let mut elements = Vec::with_capacity(total);
        let mut coefficients = Vec::with_capacity(coefficient_slots);
        enumerate_reduced_coefficients::<S>(
            coefficient_slots,
            &base_elements,
            &mut coefficients,
            &mut elements,
        );
        elements
    }
}

impl<S: ExtensionFieldSpec> SqrtField for ExtensionField<S>
where
    S::Base: EnumerableFiniteField,
{
    /// Returns one square root by exhaustive search over the full finite field.
    ///
    /// This is intentionally educational rather than asymptotically efficient.
    /// It is appropriate only for the same small finite extension backends for
    /// which [`EnumerableFiniteField`] is honest.
    fn sqrt(x: &Self::Elem) -> Option<Self::Elem> {
        Self::elements()
            .into_iter()
            .find(|candidate| Self::eq(&Self::square(candidate), x))
    }
}

impl<S: ExtensionFieldSpec> CbrtField for ExtensionField<S>
where
    S::Base: EnumerableFiniteField,
{
    /// Returns one cube root by exhaustive search over the full finite field.
    ///
    /// This is intentionally educational rather than asymptotically efficient.
    /// It is appropriate only for the same small finite extension backends for
    /// which [`EnumerableFiniteField`] is honest.
    fn cbrt(x: &Self::Elem) -> Option<Self::Elem> {
        Self::elements()
            .into_iter()
            .find(|candidate| Self::eq(&Self::cube(candidate), x))
    }
}

impl<S: ExtensionFieldSpec> QuadraticCharacterFiniteField for ExtensionField<S> where
    ExtensionField<S>: FiniteField
{
}

fn enumerate_reduced_coefficients<S: ExtensionFieldSpec>(
    remaining_slots: usize,
    base_elements: &[BaseElem<S>],
    coefficients: &mut Vec<BaseElem<S>>,
    output: &mut Vec<ExtensionFieldElement<S>>,
) where
    S::Base: EnumerableFiniteField,
{
    if remaining_slots == 0 {
        output.push(ExtensionFieldElement::<S>::new(coefficients.clone()));
        return;
    }

    for coefficient in base_elements {
        coefficients.push(coefficient.clone());
        enumerate_reduced_coefficients::<S>(
            remaining_slots - 1,
            base_elements,
            coefficients,
            output,
        );
        coefficients.pop();
    }
}
