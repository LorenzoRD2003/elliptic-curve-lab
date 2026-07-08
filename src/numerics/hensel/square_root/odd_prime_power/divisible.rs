use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};

use crate::numerics::{hensel::HenselLiftError, positive_mod_biguint};

use super::{unit::sqrt_mod_odd_prime_power, validation::validate_divisible_input};

/// Returns all square roots of `value` modulo `p^e` when `p | value`.
///
/// This is the complementary odd-prime route to `sqrt_mod_odd_prime_power`.
/// It handles the singular case by separating the `p`-adic valuation of the
/// radicand:
///
/// - if `value ≡ 0 mod p^e`, roots are exactly the multiples of
///   `p^ceil(e / 2)` modulo `p^e`
/// - if `v_p(value) = 2t < e`, solve the unit equation after dividing by
///   `p^(2t)`, then expand each unit root into `p^t` roots modulo `p^e`
/// - if `v_p(value)` is odd, no square root exists
///
/// The output contains sorted canonical representatives modulo `p^e`.
pub(crate) fn sqrt_mod_odd_prime_power_divisible_radicand(
    value: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<Vec<BigUint>, HenselLiftError> {
    validate_divisible_input(value, p, e)?;

    let modulus = p.pow(e);
    let value_modulus = positive_mod_biguint(value, &modulus);
    let valuation = p_adic_valuation_capped(&value_modulus, p, e);

    if valuation == e {
        return Ok(roots_of_zero_mod_prime_power(p, e));
    }
    if valuation % 2 == 1 {
        return Err(HenselLiftError::NoSquareRootModuloPrimePower);
    }

    let scale_exponent = valuation / 2;
    let reduced_exponent = e - valuation;
    let reduced_value = BigInt::from(&value_modulus / p.pow(valuation));
    let unit_roots = sqrt_mod_odd_prime_power(&reduced_value, p, reduced_exponent).map_err(
        |error| match error {
            HenselLiftError::QuadraticNonResidueModPrime => {
                HenselLiftError::NoSquareRootModuloPrimePower
            }
            other => other,
        },
    )?;
    let roots = expand_divisible_radicand_roots(
        [unit_roots.0, unit_roots.1],
        p,
        e,
        scale_exponent,
        reduced_exponent,
    );
    Ok(roots)
}

fn p_adic_valuation_capped(value: &BigUint, p: &BigUint, cap: u32) -> u32 {
    if value.is_zero() {
        return cap;
    }

    let mut valuation = 0u32;
    let mut residual = value.clone();
    while valuation < cap && (&residual % p).is_zero() {
        residual /= p;
        valuation += 1;
    }
    valuation
}

fn roots_of_zero_mod_prime_power(p: &BigUint, e: u32) -> Vec<BigUint> {
    let modulus = p.pow(e);
    let step = p.pow(e.div_ceil(2));
    let count = p.pow(e / 2);
    let mut roots = Vec::new();
    let mut multiplier = BigUint::zero();

    while multiplier < count {
        roots.push((&step * &multiplier) % &modulus);
        multiplier += BigUint::one();
    }
    roots
}

fn expand_divisible_radicand_roots(
    unit_roots: [BigUint; 2],
    p: &BigUint,
    e: u32,
    scale_exponent: u32,
    reduced_exponent: u32,
) -> Vec<BigUint> {
    let modulus = p.pow(e);
    let scale = p.pow(scale_exponent);
    let expansion_step = p.pow(e - scale_exponent);
    let expansion_count = p.pow(scale_exponent);
    let mut roots = Vec::new();

    for unit_root in unit_roots {
        let base_root = (&scale * unit_root) % &modulus;
        let mut offset_multiplier = BigUint::zero();
        while offset_multiplier < expansion_count {
            roots.push((&base_root + (&expansion_step * &offset_multiplier)) % &modulus);
            offset_multiplier += BigUint::one();
        }
    }

    debug_assert_eq!(
        expansion_step,
        &scale * p.pow(reduced_exponent),
        "expansion step should preserve the reduced unit root"
    );
    roots.sort();
    roots.dedup();
    roots
}
