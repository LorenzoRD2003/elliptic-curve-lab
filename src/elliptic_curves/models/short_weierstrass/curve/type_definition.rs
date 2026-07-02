use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::{CurveError, affine::AffinePoint};

/// Short-Weierstrass curve model `y^2 = x^3 + ax + b`.
///
/// This educational implementation currently supports only fields of
/// characteristic different from `2` and `3`, where the classical short form
/// and its discriminant formula behave as expected.
pub struct ShortWeierstrassCurve<F: Field> {
    a: F::Elem,
    b: F::Elem,
}

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Builds a validated short-Weierstrass curve descriptor.
    pub fn new(a: F::Elem, b: F::Elem) -> Result<Self, CurveError> {
        if F::has_characteristic(2) || F::has_characteristic(3) {
            let characteristic = F::characteristic().to_biguint();
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let curve = Self::from_validated_coefficients_unchecked(a, b);
        if F::is_zero(&curve.discriminant()) {
            return Err(CurveError::SingularCurve);
        }

        Ok(curve)
    }

    /// Returns the `a` coefficient in `x^3 + ax + b`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the `b` coefficient in `x^3 + ax + b`.
    pub fn b(&self) -> &F::Elem {
        &self.b
    }

    pub(super) fn from_validated_coefficients_unchecked(a: F::Elem, b: F::Elem) -> Self {
        Self { a, b }
    }

    /// Returns the short-Weierstrass equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("y^2 = x^3 + ({})x + ({})", self.a, self.b)
    }

    /// Builds a finite affine point without checking the curve equation.
    pub(crate) fn unchecked_point(&self, x: F::Elem, y: F::Elem) -> AffinePoint<F> {
        AffinePoint::new(x, y)
    }
}
