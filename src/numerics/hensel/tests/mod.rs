use num_bigint::{BigInt, BigUint};

use crate::polynomials::IntegerPolynomial;

mod integer_roots;
mod simple_lift;
mod square_root_lift;
mod square_root_modulus;
mod square_root_prime_power;

fn bi(value: i64) -> BigInt {
    BigInt::from(value)
}

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn polynomial(values: &[i64]) -> Vec<BigInt> {
    values.iter().copied().map(bi).collect()
}

fn integer_polynomial(values: &[i64]) -> IntegerPolynomial {
    IntegerPolynomial::new(polynomial(values))
}

fn brute_force_square_roots_mod_two_power(value: i64, e: u32) -> Vec<BigUint> {
    let modulus = 1u64 << e;
    brute_force_square_roots_modulus(value, modulus)
}

fn brute_force_square_roots_mod_prime_power(value: i64, p: u64, e: u32) -> Vec<BigUint> {
    let modulus = p.pow(e);
    brute_force_square_roots_modulus(value, modulus)
}

fn brute_force_square_roots_modulus(value: i64, modulus: u64) -> Vec<BigUint> {
    let modulus_u64 = modulus;
    let modulus_bigint = BigInt::from(modulus);
    let target = ((BigInt::from(value) % &modulus_bigint) + &modulus_bigint) % &modulus_bigint;
    let target = target
        .to_biguint()
        .expect("canonical residue modulo a positive modulus should be non-negative");
    let modulus = BigUint::from(modulus);

    (0..modulus_u64)
        .map(BigUint::from)
        .filter(|candidate| (candidate * candidate) % &modulus == target)
        .collect()
}
