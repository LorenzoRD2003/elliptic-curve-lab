use crate::elliptic_curves::CurveError;

/// Internal affine point shape used by the shared short-Weierstrass formulas.
///
/// This deliberately mirrors the two public point surfaces that currently use
/// the same group law:
///
/// - `AffinePoint<F>` over the base field
/// - `ShortWeierstrassFunctionFieldPoint<F>` over `F(E)`
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ShortWeierstrassFormulaPoint<T> {
    Infinity,
    Affine { x: T, y: T },
}

/// Internal coordinate operations needed by the shared affine short-Weierstrass formulas.
pub(crate) trait ShortWeierstrassFormulaCoordOps {
    type Coord: Clone;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError>;
    fn from_i64(&self, value: i64) -> Self::Coord;
    fn is_zero(&self, value: &Self::Coord) -> bool;
    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool;

    fn square(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        self.mul(value, value)
    }
}

fn divide_by_nonzero<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    numerator: &O::Coord,
    denominator: &O::Coord,
) -> Result<O::Coord, CurveError> {
    let inverse = ops.inv(denominator)?;
    ops.mul(numerator, &inverse)
}

fn slope_for_addition<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    x_left: &O::Coord,
    y_left: &O::Coord,
    x_right: &O::Coord,
    y_right: &O::Coord,
) -> Result<O::Coord, CurveError> {
    let numerator = ops.sub(y_right, y_left)?;
    let denominator = ops.sub(x_right, x_left)?;
    divide_by_nonzero(ops, &numerator, &denominator)
}

fn slope_for_doubling<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    curve_a: &O::Coord,
    x: &O::Coord,
    y: &O::Coord,
) -> Result<O::Coord, CurveError> {
    let three_x_squared = ops.mul(&ops.from_i64(3), &ops.square(x)?)?;
    let numerator = ops.add(&three_x_squared, curve_a)?;
    let denominator = ops.mul(&ops.from_i64(2), y)?;
    divide_by_nonzero(ops, &numerator, &denominator)
}

fn point_from_slope<O: ShortWeierstrassFormulaCoordOps>(
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

fn are_inverse_points<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    left: &ShortWeierstrassFormulaPoint<O::Coord>,
    right: &ShortWeierstrassFormulaPoint<O::Coord>,
) -> Result<bool, CurveError> {
    Ok(match (left, right) {
        (
            ShortWeierstrassFormulaPoint::Affine {
                x: x_left,
                y: y_left,
            },
            ShortWeierstrassFormulaPoint::Affine {
                x: x_right,
                y: y_right,
            },
        ) => ops.eq(x_left, x_right) && ops.is_zero(&ops.add(y_left, y_right)?),
        _ => false,
    })
}

pub(crate) fn add_formula_points<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    curve_a: &O::Coord,
    left: &ShortWeierstrassFormulaPoint<O::Coord>,
    right: &ShortWeierstrassFormulaPoint<O::Coord>,
) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
    match (left, right) {
        (ShortWeierstrassFormulaPoint::Infinity, _) => Ok(right.clone()),
        (_, ShortWeierstrassFormulaPoint::Infinity) => Ok(left.clone()),
        (
            ShortWeierstrassFormulaPoint::Affine {
                x: x_left,
                y: y_left,
            },
            ShortWeierstrassFormulaPoint::Affine {
                x: x_right,
                y: y_right,
            },
        ) => {
            if are_inverse_points(ops, left, right)? {
                return Ok(ShortWeierstrassFormulaPoint::Infinity);
            }

            if ops.eq(x_left, x_right) {
                return double_formula_point(ops, curve_a, left);
            }

            let slope = slope_for_addition(ops, x_left, y_left, x_right, y_right)?;
            point_from_slope(ops, &slope, x_left, y_left, x_right)
        }
    }
}

pub(crate) fn double_formula_point<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    curve_a: &O::Coord,
    point: &ShortWeierstrassFormulaPoint<O::Coord>,
) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
    match point {
        ShortWeierstrassFormulaPoint::Infinity => Ok(ShortWeierstrassFormulaPoint::Infinity),
        ShortWeierstrassFormulaPoint::Affine { x, y } => {
            if ops.is_zero(y) {
                return Ok(ShortWeierstrassFormulaPoint::Infinity);
            }

            let slope = slope_for_doubling(ops, curve_a, x, y)?;
            point_from_slope(ops, &slope, x, y, x)
        }
    }
}

pub(crate) fn mul_formula_point<O: ShortWeierstrassFormulaCoordOps>(
    ops: &O,
    curve_a: &O::Coord,
    point: &ShortWeierstrassFormulaPoint<O::Coord>,
    scalar: u64,
) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
    let mut result = ShortWeierstrassFormulaPoint::Infinity;
    let mut base = point.clone();
    let mut k = scalar;

    while k > 0 {
        if k & 1 == 1 {
            result = add_formula_points(ops, curve_a, &result, &base)?;
        }

        k >>= 1;
        if k > 0 {
            base = double_formula_point(ops, curve_a, &base)?;
        }
    }

    Ok(result)
}
