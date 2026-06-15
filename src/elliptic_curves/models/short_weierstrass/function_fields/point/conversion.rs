use crate::elliptic_curves::{
    CurveError,
    models::short_weierstrass::group_law_core::ShortWeierstrassFormulaPoint,
    short_weierstrass::function_fields::{
        ShortWeierstrassFunction, ShortWeierstrassFunctionFieldPoint,
    },
};
use crate::fields::traits::Field;

impl<F: Field> ShortWeierstrassFunctionFieldPoint<F> {
    pub(crate) fn to_formula_point(
        &self,
    ) -> ShortWeierstrassFormulaPoint<ShortWeierstrassFunction<F>> {
        match self {
            Self::Infinity => ShortWeierstrassFormulaPoint::Infinity,
            Self::Affine { x, y } => ShortWeierstrassFormulaPoint::Affine {
                x: x.clone(),
                y: y.clone(),
            },
        }
    }

    pub(crate) fn from_formula_point(
        point: ShortWeierstrassFormulaPoint<ShortWeierstrassFunction<F>>,
    ) -> Result<Self, CurveError> {
        match point {
            ShortWeierstrassFormulaPoint::Infinity => Ok(Self::Infinity),
            ShortWeierstrassFormulaPoint::Affine { x, y } => Self::affine(x, y),
        }
    }
}
