use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::models::short_weierstrass::group_law_core::{
    ShortWeierstrassFormulaOps, ShortWeierstrassFormulaPoint,
    formulas::{result_point_from_slope, slope_for_addition, slope_for_doubling},
};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Internal runner for the shared affine short-Weierstrass formulas.
///
/// The runner packages together the coordinate backend and the curve
/// coefficient `a`, since every caller always supplies them as one pair.
pub(crate) struct ShortWeierstrassFormulaRunner<'a, O: ShortWeierstrassFormulaOps> {
    ops: &'a O,
    curve_a: &'a O::Coord,
}

impl<'a, O: ShortWeierstrassFormulaOps> ShortWeierstrassFormulaRunner<'a, O> {
    pub(crate) fn new(ops: &'a O, curve_a: &'a O::Coord) -> Self {
        Self { ops, curve_a }
    }

    fn points_are_vertical_opposites(
        &self,
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
            ) => self.ops.eq(x_left, x_right) && self.ops.is_zero(&self.ops.add(y_left, y_right)?),
            _ => false,
        })
    }

    pub(crate) fn add_points(
        &self,
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
                if self.points_are_vertical_opposites(left, right)? {
                    return Ok(ShortWeierstrassFormulaPoint::Infinity);
                }

                if self.ops.eq(x_left, x_right) {
                    return self.double_point(left);
                }

                let slope = slope_for_addition(self.ops, x_left, y_left, x_right, y_right)?;
                result_point_from_slope(self.ops, &slope, x_left, y_left, x_right)
            }
        }
    }

    pub(crate) fn double_point(
        &self,
        point: &ShortWeierstrassFormulaPoint<O::Coord>,
    ) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
        match point {
            ShortWeierstrassFormulaPoint::Infinity => Ok(ShortWeierstrassFormulaPoint::Infinity),
            ShortWeierstrassFormulaPoint::Affine { x, y } => {
                if self.ops.is_zero(y) {
                    return Ok(ShortWeierstrassFormulaPoint::Infinity);
                }

                let slope = slope_for_doubling(self.ops, self.curve_a, x, y)?;
                result_point_from_slope(self.ops, &slope, x, y, x)
            }
        }
    }

    /// Multiplies one affine-formula point by a nonnegative scalar via the
    /// binary double-and-add method.
    ///
    /// Complexity: `Θ(log n)` runner additions/doublings for scalar `n`,
    /// counting only calls to the shared short-Weierstrass group-law core.
    ///
    pub(crate) fn mul_point(
        &self,
        point: &ShortWeierstrassFormulaPoint<O::Coord>,
        scalar: impl crate::elliptic_curves::traits::ScalarInput,
    ) -> Result<ShortWeierstrassFormulaPoint<O::Coord>, CurveError> {
        let mut result = ShortWeierstrassFormulaPoint::Infinity;
        let mut base = point.clone();
        let mut k = scalar.into_biguint_scalar();

        while !k.is_zero() {
            if (&k & BigUint::one()) == BigUint::one() {
                result = self.add_points(&result, &base)?;
            }

            k >>= 1usize;
            if !k.is_zero() {
                base = self.double_point(&base)?;
            }
        }

        Ok(result)
    }
}
