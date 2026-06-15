use crate::elliptic_curves::CurveError;

/// Internal coordinate operations needed by the shared affine
/// short-Weierstrass formulas.
///
/// This is intentionally smaller than the full field-family traits. It models
/// only the operations that the secant/tangent formulas actually use, and lets
/// the same core arithmetic run both over ordinary field coordinates and over
/// the function-field value layer.
pub(crate) trait ShortWeierstrassFormulaOps {
    type Coord: Clone;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn lift_i64(&self, value: i64) -> Self::Coord;
    fn is_zero(&self, value: &Self::Coord) -> bool;
    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool;

    fn div(
        &self,
        numerator: &Self::Coord,
        denominator: &Self::Coord,
    ) -> Result<Self::Coord, CurveError> {
        let inverse = self.inv(denominator)?;
        self.mul(numerator, &inverse)
    }

    fn square(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        self.mul(value, value)
    }
}
