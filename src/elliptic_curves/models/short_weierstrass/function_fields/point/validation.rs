use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    short_weierstrass::function_fields::{
        ShortWeierstrassFunction, ShortWeierstrassFunctionFieldPoint,
    },
};
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassFunctionFieldPoint<F> {
    pub(crate) fn validate_affine_coordinates_on_curve(
        x: &ShortWeierstrassFunction<F>,
        y: &ShortWeierstrassFunction<F>,
        curve: &ShortWeierstrassCurve<F>,
    ) -> Result<(), CurveError> {
        if !F::eq(x.curve().a(), curve.a()) || !F::eq(x.curve().b(), curve.b()) {
            return Err(CurveError::IncompatibleFunctionFieldCurves);
        }
        if !F::eq(y.curve().a(), curve.a()) || !F::eq(y.curve().b(), curve.b()) {
            return Err(CurveError::IncompatibleFunctionFieldCurves);
        }

        let lhs = y.mul(y)?;
        let rhs = curve.evaluate_curve_rhs_at_function_x(x)?;
        if lhs == rhs {
            Ok(())
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }

    pub(crate) fn ensure_point_lives_on(
        &self,
        curve: &ShortWeierstrassCurve<F>,
    ) -> Result<(), CurveError> {
        match self {
            Self::Infinity => Ok(()),
            Self::Affine { x, y } => Self::validate_affine_coordinates_on_curve(x, y, curve),
        }
    }
}
