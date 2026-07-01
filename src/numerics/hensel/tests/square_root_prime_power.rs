use crate::numerics::hensel::{
    HenselLiftError, sqrt_mod_odd_prime_power, sqrt_mod_odd_prime_power_divisible_radicand,
    sqrt_mod_two_power,
};

use super::{
    bi, brute_force_square_roots_mod_prime_power, brute_force_square_roots_mod_two_power, bu,
};

#[test]
fn sqrt_mod_odd_prime_power_returns_both_canonical_roots() {
    let roots = sqrt_mod_odd_prime_power(&bi(2), &bu(7), 4)
        .expect("2 should be a square modulo 7 and lift to 7^4");
    let modulus = bu(7).pow(4);

    assert_eq!(roots, (bu(235), bu(2166)));
    assert_eq!((&roots.0 * &roots.0 - bu(2)) % &modulus, bu(0));
    assert_eq!((&roots.1 * &roots.1 - bu(2)) % &modulus, bu(0));
    assert_eq!((&roots.0 + &roots.1) % &modulus, bu(0));
}

#[test]
fn sqrt_mod_odd_prime_power_normalizes_negative_radicands() {
    let roots =
        sqrt_mod_odd_prime_power(&bi(-1), &bu(5), 3).expect("-1 should be a square modulo 5^3");
    let modulus = bu(5).pow(3);

    assert_eq!(roots, (bu(57), bu(68)));
    assert_eq!((&roots.0 * &roots.0 + bu(1)) % &modulus, bu(0));
    assert_eq!((&roots.1 * &roots.1 + bu(1)) % &modulus, bu(0));
}

#[test]
fn sqrt_mod_odd_prime_power_rejects_postponed_or_impossible_cases() {
    assert_eq!(
        sqrt_mod_odd_prime_power(&bi(2), &bu(7), 0),
        Err(HenselLiftError::ZeroTargetLevel)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power(&bi(1), &bu(2), 3),
        Err(HenselLiftError::EvenPrimeUnsupported)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power(&bi(7), &bu(7), 2),
        Err(HenselLiftError::RadicandDivisibleByPrimeUnsupported)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power(&bi(3), &bu(7), 2),
        Err(HenselLiftError::QuadraticNonResidueModPrime)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power(&bi(1), &bu(9), 2),
        Err(HenselLiftError::NonPrimeModulus)
    );
}

#[test]
fn sqrt_mod_odd_prime_power_divisible_handles_zero_radicand() {
    let roots = sqrt_mod_odd_prime_power_divisible_radicand(&bi(0), &bu(3), 4)
        .expect("0 should have roots modulo 3^4");

    assert_eq!(
        roots,
        vec![
            bu(0),
            bu(9),
            bu(18),
            bu(27),
            bu(36),
            bu(45),
            bu(54),
            bu(63),
            bu(72)
        ]
    );
}

#[test]
fn sqrt_mod_odd_prime_power_divisible_handles_even_valuation() {
    let roots = sqrt_mod_odd_prime_power_divisible_radicand(&bi(9), &bu(3), 4)
        .expect("9 should have roots modulo 3^4");
    let modulus = bu(3).pow(4);

    assert_eq!(roots, vec![bu(3), bu(24), bu(30), bu(51), bu(57), bu(78)]);
    for root in roots {
        assert_eq!((&root * &root - bu(9)) % &modulus, bu(0));
    }
}

#[test]
fn sqrt_mod_odd_prime_power_divisible_rejects_odd_valuation_and_wrong_route() {
    assert_eq!(
        sqrt_mod_odd_prime_power_divisible_radicand(&bi(3), &bu(3), 4),
        Err(HenselLiftError::NoSquareRootModuloPrimePower)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power_divisible_radicand(&bi(1), &bu(3), 4),
        Err(HenselLiftError::RadicandNotDivisibleByPrime)
    );
    assert_eq!(
        sqrt_mod_odd_prime_power_divisible_radicand(&bi(0), &bu(2), 4),
        Err(HenselLiftError::EvenPrimeUnsupported)
    );
}

#[test]
fn sqrt_mod_odd_prime_power_divisible_matches_brute_force_for_small_cases() {
    for e in 1..=5 {
        for value in [-27, -18, -9, -6, -3, 0, 3, 6, 9, 18, 27] {
            let expected = brute_force_square_roots_mod_prime_power(value, 3, e);
            let actual = sqrt_mod_odd_prime_power_divisible_radicand(&bi(value), &bu(3), e);

            if expected.is_empty() {
                assert_eq!(actual, Err(HenselLiftError::NoSquareRootModuloPrimePower));
            } else {
                assert_eq!(actual.expect("brute-force roots should exist"), expected);
            }
        }
    }
}

#[test]
fn sqrt_mod_two_power_handles_small_structural_cases() {
    assert_eq!(
        sqrt_mod_two_power(&bi(0), 1).expect("0 has a square root modulo 2"),
        vec![bu(0)]
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(1), 1).expect("1 has a square root modulo 2"),
        vec![bu(1)]
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(0), 2).expect("0 has roots modulo 4"),
        vec![bu(0), bu(2)]
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(1), 2).expect("1 has roots modulo 4"),
        vec![bu(1), bu(3)]
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(1), 3).expect("1 has roots modulo 8"),
        vec![bu(1), bu(3), bu(5), bu(7)]
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(4), 4).expect("4 has roots modulo 16"),
        vec![bu(2), bu(6), bu(10), bu(14)]
    );
}

#[test]
fn sqrt_mod_two_power_rejects_zero_exponent_and_non_squares() {
    assert_eq!(
        sqrt_mod_two_power(&bi(1), 0),
        Err(HenselLiftError::ZeroTargetLevel)
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(2), 2),
        Err(HenselLiftError::NoSquareRootModuloPrimePower)
    );
    assert_eq!(
        sqrt_mod_two_power(&bi(-1), 3),
        Err(HenselLiftError::NoSquareRootModuloPrimePower)
    );
}

#[test]
fn sqrt_mod_two_power_matches_brute_force_for_small_exponents() {
    for e in 1..=7 {
        for value in -4..=20 {
            let expected = brute_force_square_roots_mod_two_power(value, e);
            let actual = sqrt_mod_two_power(&bi(value), e);

            if expected.is_empty() {
                assert_eq!(actual, Err(HenselLiftError::NoSquareRootModuloPrimePower));
            } else {
                assert_eq!(actual.expect("brute-force roots should lift"), expected);
            }
        }
    }
}
