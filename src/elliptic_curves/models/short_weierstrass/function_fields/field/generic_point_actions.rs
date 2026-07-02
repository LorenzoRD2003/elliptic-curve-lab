use crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField;
use crate::elliptic_curves::{
    AffinePoint, CurveError, short_weierstrass::function_fields::ShortWeierstrassFunctionFieldPoint,
};
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassFunctionField<F> {
    pub(crate) fn translate_generic_point_by_base_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        let generic_point = self.generic_point();
        let constant_point = self.embed_constant_affine_point(point)?;
        self.add_points(&generic_point, &constant_point)
    }

    #[cfg(test)]
    pub(crate) fn double_generic_point(
        &self,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.double_point(&self.generic_point())
    }

    pub(crate) fn generic_point_multiple(
        &self,
        scalar: u64,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.mul_scalar_point(&self.generic_point(), scalar)
    }
}
