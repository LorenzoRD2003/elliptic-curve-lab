use core::marker::PhantomData;

use crate::fields::{FieldError, rational_function_field::RationalFunction, traits::Field};
use crate::polynomials::DensePolynomial;

/// Zero-sized field family for the univariate rational function field `F(x)`.
///
/// This type separates the field-family metadata and trait implementation from
/// the stored rational-function values, which live in [`RationalFunction`].
pub struct RationalFunctionField<F: Field> {
    _base: PhantomData<F>,
}

impl<F: Field> RationalFunctionField<F> {
    /// Embeds a polynomial into the rational function field.
    pub fn from_polynomial(polynomial: DensePolynomial<F>) -> RationalFunction<F> {
        RationalFunction::<F>::from_polynomial(polynomial)
    }

    /// Builds a constant rational function.
    pub fn constant(value: F::Elem) -> RationalFunction<F> {
        RationalFunction::<F>::constant(value)
    }

    /// Returns the distinguished indeterminate `x`.
    pub fn indeterminate() -> RationalFunction<F> {
        RationalFunction::<F>::indeterminate()
    }
}

impl<F: Field> Field for RationalFunctionField<F> {
    const IS_ALGEBRAICALLY_CLOSED: bool = false;

    type Elem = RationalFunction<F>;

    fn characteristic() -> u64 {
        F::characteristic()
    }

    fn zero() -> Self::Elem {
        Self::constant(F::zero())
    }

    fn one() -> Self::Elem {
        Self::constant(F::one())
    }

    fn from_i64(n: i64) -> Self::Elem {
        Self::constant(F::from_i64(n))
    }

    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x.add(y)
    }

    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x.sub(y)
    }

    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        x.mul(y)
    }

    fn neg(x: &Self::Elem) -> Self::Elem {
        x.neg()
    }

    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        x.inverse().ok()
    }

    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }

    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        x.inverse()
    }

    fn elem_from_u64(value: u64) -> Self::Elem {
        Self::constant(F::elem_from_u64(value))
    }
}
