use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassFunction};
use crate::fields::{Field, RationalFunction};
use crate::polynomials::DensePolynomial;

pub(super) fn pow_rational_function_u128<F: Field>(
    function: &RationalFunction<F>,
    exponent: u128,
) -> RationalFunction<F> {
    let mut result = RationalFunction::<F>::constant(F::one());
    let mut base = function.clone();
    let mut exp = exponent;

    while exp > 0 {
        if exp & 1 == 1 {
            result = result.mul(&base);
        }

        exp >>= 1;

        if exp > 0 {
            base = base.mul(&base);
        }
    }

    result
}

pub(super) fn x_pullback_from_power<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    power: u128,
) -> ShortWeierstrassFunction<F> {
    ShortWeierstrassFunction::<F>::from_rational_function(
        curve.clone(),
        pow_rational_function_u128(&RationalFunction::<F>::indeterminate(), power),
    )
}

pub(super) fn y_pullback_from_power<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    power: u128,
) -> ShortWeierstrassFunction<F> {
    let rhs = RationalFunction::<F>::from_polynomial(DensePolynomial::<F>::new(vec![
        curve.b().clone(),
        curve.a().clone(),
        F::zero(),
        F::one(),
    ]));

    ShortWeierstrassFunction::<F>::new(
        curve.clone(),
        RationalFunction::<F>::constant(F::zero()),
        pow_rational_function_u128(&rhs, (power - 1) / 2),
    )
}
