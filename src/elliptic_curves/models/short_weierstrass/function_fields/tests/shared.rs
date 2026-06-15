use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{Fp, Q, rational_function_field::RationalFunction, traits::Field};
use crate::polynomials::DensePolynomial;

pub(super) type F17 = Fp<17>;

pub(super) fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
    DensePolynomial::<F17>::new(values.iter().copied().map(F17::elem_from_u64).collect())
}

pub(super) fn q_dense(values: &[(i64, i64)]) -> DensePolynomial<Q> {
    DensePolynomial::<Q>::new(
        values
            .iter()
            .map(|&(numerator, denominator)| {
                Q::div(&Q::from_i64(numerator), &Q::from_i64(denominator))
                    .expect("denominator should be non-zero")
            })
            .collect(),
    )
}

pub(super) fn f17_curve() -> ShortWeierstrassCurve<F17> {
    ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
        .expect("curve should be nonsingular")
}

pub(super) fn q_curve() -> ShortWeierstrassCurve<Q> {
    ShortWeierstrassCurve::<Q>::new(Q::from_i64(-1), Q::from_i64(0)).expect("curve should exist")
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
