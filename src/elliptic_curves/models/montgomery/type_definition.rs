use core::fmt;

use crate::elliptic_curves::CurveError;
use crate::fields::traits::Field;

/// Montgomery curve model `B y^2 = x^3 + A x^2 + x`.
pub struct MontgomeryCurve<F: Field> {
    a: F::Elem,
    b: F::Elem,
}

impl<F: Field> MontgomeryCurve<F> {
    /// Builds a validated Montgomery curve descriptor.
    pub fn new(a: F::Elem, b: F::Elem) -> Result<Self, CurveError> {
        let characteristic = F::characteristic();
        if characteristic == 2 {
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }
        if F::is_zero(&b) {
            return Err(CurveError::SingularCurve);
        }

        let curve = Self::from_validated_coefficients_unchecked(a, b);
        if F::is_zero(&curve.discriminant()) {
            return Err(CurveError::SingularCurve);
        }

        Ok(curve)
    }

    /// Returns the coefficient `A`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the coefficient `B`.
    pub fn b(&self) -> &F::Elem {
        &self.b
    }

    pub(super) fn from_validated_coefficients_unchecked(a: F::Elem, b: F::Elem) -> Self {
        Self { a, b }
    }

    /// Returns the defining equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("({})y^2 = x^3 + ({})x^2 + x", self.b, self.a)
    }
}
