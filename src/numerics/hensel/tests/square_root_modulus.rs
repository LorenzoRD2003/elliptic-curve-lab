use crate::numerics::hensel::{HenselLiftError, sqrt_mod_m};

use super::{bi, brute_force_square_roots_modulus, bu};

#[test]
fn sqrt_mod_general_modulus_combines_prime_power_roots_with_crt() {
    let roots = sqrt_mod_m(&bi(1), &bu(72)).expect("1 should have roots modulo 72 = 2^3 * 3^2");

    assert_eq!(roots, brute_force_square_roots_modulus(1, 72));
    assert_eq!(roots.len(), 8);
}

#[test]
fn sqrt_mod_general_modulus_handles_local_divisible_radicands() {
    let roots = sqrt_mod_m(&bi(9), &bu(45)).expect("9 should have roots modulo 45 = 3^2 * 5");

    assert_eq!(roots, brute_force_square_roots_modulus(9, 45));
    assert_eq!(roots.len(), 6);
}

#[test]
fn sqrt_mod_general_modulus_rejects_trivial_moduli_and_local_non_squares() {
    assert_eq!(
        sqrt_mod_m(&bi(1), &bu(1)),
        Err(HenselLiftError::TrivialModulus)
    );
    assert_eq!(
        sqrt_mod_m(&bi(5), &bu(21)),
        Err(HenselLiftError::NoSquareRootModuloPrimePower)
    );
}
