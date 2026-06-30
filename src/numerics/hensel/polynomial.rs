use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

pub(super) fn evaluate_polynomial(coefs: &[BigInt], x: &BigInt) -> BigInt {
    coefs
        .iter()
        .rev()
        .fold(BigInt::zero(), |acc, c| acc * x + c)
}

pub(super) fn evaluate_derivative(coefs: &[BigInt], x: &BigInt) -> BigInt {
    coefs
        .iter()
        .enumerate()
        .skip(1)
        .rev()
        .fold(BigInt::zero(), |acc, (deg, c)| {
            acc * x + c * BigInt::from(deg)
        })
}

pub(super) fn positive_mod_biguint(value: &BigInt, modulus: &BigUint) -> BigUint {
    positive_mod_bigint(value, &BigInt::from(modulus.clone()))
        .to_biguint()
        .expect("positive modular representative should convert to BigUint")
}

pub(super) fn positive_mod_bigint(value: &BigInt, modulus: &BigInt) -> BigInt {
    debug_assert!(modulus > &BigInt::zero());
    ((value % modulus) + modulus) % modulus
}
