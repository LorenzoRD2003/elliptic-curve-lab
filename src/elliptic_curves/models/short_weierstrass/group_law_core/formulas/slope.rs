use crate::elliptic_curves::{
    CurveError, models::short_weierstrass::group_law_core::ShortWeierstrassFormulaOps,
};

pub(crate) fn slope_for_addition<O: ShortWeierstrassFormulaOps>(
    ops: &O,
    x_left: &O::Coord,
    y_left: &O::Coord,
    x_right: &O::Coord,
    y_right: &O::Coord,
) -> Result<O::Coord, CurveError> {
    let numerator = ops.sub(y_right, y_left)?;
    let denominator = ops.sub(x_right, x_left)?;
    ops.div(&numerator, &denominator)
}

pub(crate) fn slope_for_doubling<O: ShortWeierstrassFormulaOps>(
    ops: &O,
    curve_a: &O::Coord,
    x: &O::Coord,
    y: &O::Coord,
) -> Result<O::Coord, CurveError> {
    let three_x_squared = ops.mul(&ops.lift_i64(3), &ops.square(x)?)?;
    let numerator = ops.add(&three_x_squared, curve_a)?;
    let denominator = ops.mul(&ops.lift_i64(2), y)?;
    ops.div(&numerator, &denominator)
}
