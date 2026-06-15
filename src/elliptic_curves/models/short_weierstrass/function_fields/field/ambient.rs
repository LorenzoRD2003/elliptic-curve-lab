use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{
    rational_function_field::RationalFunction,
    traits::{AmbientField, Field},
};

/// Runtime ambient family for the function field of one concrete
/// short-Weierstrass curve.
///
/// Unlike the base-field trait [`Field`], this ambient object
/// depends on a specific curve chosen at runtime.
#[derive(Clone)]
pub struct ShortWeierstrassFunctionField<F: Field> {
    curve: ShortWeierstrassCurve<F>,
}

impl<F: Field> ShortWeierstrassFunctionField<F> {
    /// Builds the ambient function field of a short Weierstrass curve.
    pub fn new(curve: ShortWeierstrassCurve<F>) -> Self {
        Self { curve }
    }

    /// Returns the ambient short-Weierstrass curve.
    pub fn curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.curve
    }

    /// Returns the zero function.
    pub fn zero(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::zero(self.curve.clone())
    }

    /// Returns the constant function `1`.
    pub fn one(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::one(self.curve.clone())
    }

    /// Returns the coordinate function `x`.
    pub fn x(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::from_rational_function(
            self.curve.clone(),
            RationalFunction::<F>::indeterminate(),
        )
    }

    /// Returns the coordinate function `y`.
    pub fn y(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::y(self.curve.clone())
    }
}

impl<F: Field> AmbientField for ShortWeierstrassFunctionField<F> {
    type Elem = ShortWeierstrassFunction<F>;
    type Error = CurveError;

    fn zero(&self) -> Self::Elem {
        Self::zero(self)
    }

    fn one(&self) -> Self::Elem {
        Self::one(self)
    }

    fn eq(&self, left: &Self::Elem, right: &Self::Elem) -> bool {
        left == right
    }

    fn add(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        left.add(right)
    }

    fn neg(&self, value: &Self::Elem) -> Self::Elem {
        value.neg()
    }

    fn mul(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        left.mul(right)
    }

    fn inverse(&self, value: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        value.inverse()
    }
}
