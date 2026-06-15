use crate::elliptic_curves::{
    CurveError,
    models::short_weierstrass::group_law_core::{
        ShortWeierstrassFormulaOps, ShortWeierstrassFormulaPoint,
    },
};

pub(crate) fn result_point_from_slope<O: ShortWeierstrassFormulaOps>(
    ops: &O,
    slope: &O::Coord,
    x_left: &O::Coord,
    y_left: &O::Coord,
    x_right: &O::Coord,
) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
    let x_result = ops.sub(&ops.sub(&ops.square(slope)?, x_left)?, x_right)?;
    let y_result = ops.sub(&ops.mul(slope, &ops.sub(x_left, &x_result)?)?, y_left)?;
    Ok(ShortWeierstrassFormulaPoint::Affine {
        x: x_result,
        y: y_result,
    })
}
