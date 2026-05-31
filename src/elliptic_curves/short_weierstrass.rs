use crate::elliptic_curves::{AffinePoint, traits::CurveEquation};

/// Tentative short-Weierstrass model `y^2 = x^3 + ax + b`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShortWeierstrassCurve<Scalar> {
    pub a: Scalar,
    pub b: Scalar,
}

impl<Scalar> ShortWeierstrassCurve<Scalar> {
    /// Builds a curve descriptor from coefficients.
    pub fn new(a: Scalar, b: Scalar) -> Self {
        Self { a, b }
    }
}

impl<Scalar> CurveEquation<AffinePoint<Scalar>> for ShortWeierstrassCurve<Scalar> {
    fn is_on_curve(&self, _point: &AffinePoint<Scalar>) -> bool {
        todo!("curve membership checks will be implemented once field arithmetic is available")
    }
}
