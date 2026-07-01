use num_bigint::{BigInt, BigUint};

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

fn brute_force_square_roots_mod_two_power(value: i64, e: u32) -> Vec<BigUint> {
    let modulus = 1u64 << e;
    brute_force_square_roots_modulus(value, modulus)
}

fn brute_force_square_roots_mod_prime_power(value: i64, p: u64, e: u32) -> Vec<BigUint> {
    let modulus = p.pow(e);
    brute_force_square_roots_modulus(value, modulus)
}

fn brute_force_square_roots_modulus(value: i64, modulus: u64) -> Vec<BigUint> {
    let target = ((value % modulus as i64) + modulus as i64) as u64 % modulus;

    (0..modulus)
        .filter(|candidate| (candidate * candidate) % modulus == target)
        .map(bu)
        .collect()
}
