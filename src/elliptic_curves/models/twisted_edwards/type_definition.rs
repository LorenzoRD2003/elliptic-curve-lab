use core::fmt;

use crate::elliptic_curves::CurveError;
use crate::fields::traits::Field;

/// Twisted Edwards curve model `E_{a,d} = a x^2 + y^2 = 1 + d x^2 y^2`.
pub struct TwistedEdwardsCurve<F: Field> {
    a: F::Elem,
    d: F::Elem,
}

impl<F: Field> TwistedEdwardsCurve<F> {
    /// Builds a validated twisted-Edwards curve descriptor.
    pub fn new(a: F::Elem, d: F::Elem) -> Result<Self, CurveError> {
        let characteristic = F::characteristic();
        if characteristic == 2 {
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let curve = Self::from_validated_coefficients_unchecked(a, d);
        if F::is_zero(&curve.discriminant()) {
            return Err(CurveError::SingularCurve);
        }

        Ok(curve)
    }

    /// Returns the coefficient `a`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the coefficient `d`.
    pub fn d(&self) -> &F::Elem {
        &self.d
    }

    pub(super) fn from_validated_coefficients_unchecked(a: F::Elem, d: F::Elem) -> Self {
        Self { a, d }
    }

    /// Returns the defining equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("({})x^2 + y^2 = 1 + ({})x^2y^2", self.a, self.d)
    }
}
