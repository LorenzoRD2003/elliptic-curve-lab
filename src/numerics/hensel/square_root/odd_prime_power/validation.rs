use num_bigint::{BigInt, BigUint};
use num_prime::nt_funcs::is_prime;
use num_traits::Zero;

use crate::numerics::{hensel::HenselLiftError, positive_mod_biguint};

pub(super) fn validate_unit_input(
    value: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<(), HenselLiftError> {
    validate_odd_prime_power(p, e)?;
    if positive_mod_biguint(value, p).is_zero() {
        return Err(HenselLiftError::RadicandDivisibleByPrimeUnsupported);
    }
    Ok(())
}

pub(super) fn validate_divisible_input(
    value: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<(), HenselLiftError> {
    validate_odd_prime_power(p, e)?;
    if !positive_mod_biguint(value, p).is_zero() {
        return Err(HenselLiftError::RadicandNotDivisibleByPrime);
    }
    Ok(())
}

fn validate_odd_prime_power(p: &BigUint, e: u32) -> Result<(), HenselLiftError> {
    if e == 0 {
        return Err(HenselLiftError::ZeroTargetLevel);
    }
    if p < &BigUint::from(2u8) || !is_prime(p, None).probably() {
        return Err(HenselLiftError::NonPrimeModulus);
    }
    if p == &BigUint::from(2u8) {
        return Err(HenselLiftError::EvenPrimeUnsupported);
    }
    Ok(())
}
