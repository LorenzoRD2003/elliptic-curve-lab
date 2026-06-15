use crate::elliptic_curves::{AffinePoint, CurveError, ShortWeierstrassCurve, traits::CurveModel};
use crate::fields::traits::FiniteField;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Applies the relative Frobenius `π_q` to a point on this curve over the
    /// represented base field `F_q`.
    ///
    /// In the current backend model this action is the identity on
    /// coordinates, so this returns the validated input point.
    ///
    /// Complexity: Θ(1)
    pub(crate) fn relative_frobenius_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        Ok(point.clone())
    }
}
