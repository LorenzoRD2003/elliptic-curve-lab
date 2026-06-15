use crate::elliptic_curves::{
    AffinePoint, CurveError,
    short_weierstrass::function_fields::{
        ShortWeierstrassFunction, ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
    },
    traits::CurveModel,
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassFunctionField<F> {
    /// Embeds one affine base-field point as a constant point of `F(E)`.
    pub(crate) fn embed_constant_affine_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        if !self.curve().contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(ShortWeierstrassFunctionFieldPoint::Infinity),
            AffinePoint::Finite { .. } => ShortWeierstrassFunctionFieldPoint::affine(
                ShortWeierstrassFunction::from_finite_affine_point(
                    self.curve().clone(),
                    point,
                    true,
                ),
                ShortWeierstrassFunction::from_finite_affine_point(
                    self.curve().clone(),
                    point,
                    false,
                ),
            ),
        }
    }
}
