use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::numerics::perfect_powers::PerfectPowerOutcome;

pub(super) fn validate_search_input(input: &BigUint) -> Result<u32, PerfectPowerOutcome> {
    if input <= &BigUint::one() {
        return Err(PerfectPowerOutcome::DegenerateInput);
    }
    if !is_coprime_to_six(input) {
        return Err(PerfectPowerOutcome::NotCoprimeToSix);
    }
    max_prime_exponent(input).ok_or(PerfectPowerOutcome::ExponentTooLarge)
}

fn is_coprime_to_six(input: &BigUint) -> bool {
    !(input % BigUint::from(2u8)).is_zero() && !(input % BigUint::from(3u8)).is_zero()
}

fn max_prime_exponent(input: &BigUint) -> Option<u32> {
    let floor_log2 = input.bits().checked_sub(1)?;
    u32::try_from(floor_log2).ok()
}
