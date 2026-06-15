use proptest::prelude::*;

use crate::fields::{Fp, traits::Field};
use crate::{
    elliptic_curves::{
        AffinePoint, ShortWeierstrassCurve,
        short_weierstrass::function_fields::ShortWeierstrassFunction,
        traits::{CurveModel, EnumerableCurveModel},
    },
    fields::rational_function_field::RationalFunction,
};

pub(super) type F41 = Fp<41>;
pub(super) type Curve = ShortWeierstrassCurve<F41>;

pub(super) fn curve() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

pub(super) fn curve_and_point() -> impl Strategy<Value = (Curve, <Curve as CurveModel>::Point)> {
    let curve = curve();
    let points = curve.points();
    let len = points.len();

    (0usize..len).prop_map(move |index| (curve.clone(), points[index].clone()))
}

pub(super) fn evaluate_short_weierstrass_function_at_point<F: Field>(
    function: &ShortWeierstrassFunction<F>,
    point: &AffinePoint<F>,
) -> Option<F::Elem> {
    let x = AffinePoint::x_coordinate(point)?;
    let y = AffinePoint::y_coordinate(point)?;
    let a_value = evaluate_rational_function_at_x(function.a_part(), x)?;
    let b_value = evaluate_rational_function_at_x(function.b_part(), x)?;

    Some(F::add(&a_value, &F::mul(y, &b_value)))
}

fn evaluate_rational_function_at_x<F: Field>(
    function: &RationalFunction<F>,
    x: &F::Elem,
) -> Option<F::Elem> {
    let numerator = function.numerator().evaluate(x).ok()?;
    let denominator = function.denominator().evaluate(x).ok()?;

    if F::is_zero(&denominator) {
        None
    } else {
        F::div(&numerator, &denominator).ok()
    }
}
