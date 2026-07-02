use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::rational_function_field::RationalFunction;
use crate::fields::traits::*;
use crate::polynomials::DensePolynomial;
use num_bigint::BigUint;

pub(super) fn x_pullback_from_power<F: Field>(
    curve: &ShortWeierstrassCurve<F>,
    power: u128,
) -> ShortWeierstrassFunction<F> {
    ShortWeierstrassFunction::<F>::from_rational_function(
        curve.clone(),
        RationalFunction::<F>::indeterminate().pow_biguint(&BigUint::from(power)),
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
        rhs.pow_biguint(&((BigUint::from(power) - BigUint::from(1u8)) / BigUint::from(2u8))),
    )
}
