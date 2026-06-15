use crate::elliptic_curves::{
    CurveError,
    models::short_weierstrass::group_law_core::{
        ShortWeierstrassFormulaOps, ShortWeierstrassFormulaPoint, ShortWeierstrassFormulaRunner,
    },
};
use crate::fields::{Q, traits::Field};

struct RationalOps;

impl ShortWeierstrassFormulaOps for RationalOps {
    type Coord = <Q as Field>::Elem;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(Q::add(left, right))
    }

    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(Q::sub(left, right))
    }

    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(Q::mul(left, right))
    }

    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Q::inverse(value).map_err(|_| CurveError::NonInvertibleFunctionFieldElement)
    }

    fn lift_i64(&self, value: i64) -> Self::Coord {
        Q::from_i64(value)
    }

    fn is_zero(&self, value: &Self::Coord) -> bool {
        Q::eq(value, &Q::zero())
    }

    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool {
        Q::eq(left, right)
    }
}

fn affine(x: i64, y: i64) -> ShortWeierstrassFormulaPoint<<Q as Field>::Elem> {
    ShortWeierstrassFormulaPoint::Affine {
        x: Q::from_i64(x),
        y: Q::from_i64(y),
    }
}

fn curve_a() -> <Q as Field>::Elem {
    Q::from_i64(-1)
}

fn runner() -> ShortWeierstrassFormulaRunner<'static, RationalOps> {
    static OPS: RationalOps = RationalOps;
    let curve_a = Box::leak(Box::new(curve_a()));
    ShortWeierstrassFormulaRunner::new(&OPS, curve_a)
}

#[test]
fn add_points_respects_the_identity() {
    let runner = runner();
    let point = affine(0, 0);

    assert_eq!(
        runner.add_points(&ShortWeierstrassFormulaPoint::Infinity, &point),
        Ok(point.clone())
    );
    assert_eq!(
        runner.add_points(&point, &ShortWeierstrassFormulaPoint::Infinity),
        Ok(point)
    );
}

#[test]
fn add_points_detects_vertical_opposites() {
    let runner = runner();
    let point = affine(0, 0);

    assert_eq!(
        runner.add_points(&point, &point),
        Ok(ShortWeierstrassFormulaPoint::Infinity)
    );
}

#[test]
fn double_point_returns_infinity_when_y_is_zero() {
    let runner = runner();

    assert_eq!(
        runner.double_point(&affine(0, 0)),
        Ok(ShortWeierstrassFormulaPoint::Infinity)
    );
}

#[test]
fn doubling_matches_adding_a_point_to_itself() {
    let runner = runner();
    let point = affine(1, 1);

    assert_eq!(
        runner.double_point(&point),
        runner.add_points(&point, &point)
    );
}

#[test]
fn scalar_multiplication_handles_small_scalars() {
    let runner = runner();
    let point = affine(1, 1);

    assert_eq!(
        runner.mul_point(&point, 0),
        Ok(ShortWeierstrassFormulaPoint::Infinity)
    );
    assert_eq!(runner.mul_point(&point, 1), Ok(point.clone()));
    assert_eq!(runner.mul_point(&point, 2), runner.double_point(&point));
    assert_eq!(
        runner.mul_point(&point, 3),
        runner.add_points(
            &runner.double_point(&point).expect("2P should exist"),
            &point
        )
    );
}
