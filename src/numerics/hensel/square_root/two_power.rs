use num_bigint::{BigInt, BigUint};

use crate::numerics::hensel::{HenselLiftError, polynomial::positive_mod_biguint};

/// Returns all square roots of `value` modulo `2^e`.
///
/// Unlike the odd-prime simple-root route, the derivative `2x` is never a unit
/// modulo `2`. This helper therefore uses direct bit-by-bit lifting: each root
/// modulo `2^k` has two possible representatives modulo `2^(k + 1)`, and the
/// function keeps exactly those whose square matches `value`.
///
/// The output contains sorted canonical representatives modulo `2^e`. If no
/// lift survives, the function returns [`HenselLiftError::NoSquareRootModuloPrimePower`].
///
/// Complexity: `Θ(e * R)` modular squarings, where `R` is the total number of
/// candidate roots inspected across all lifting levels. The output itself may
/// contain many roots for highly divisible radicands, so this route is intended
/// as exact shared infrastructure rather than a constant-size helper.
pub(crate) fn sqrt_mod_two_power(value: &BigInt, e: u32) -> Result<Vec<BigUint>, HenselLiftError> {
    if e == 0 {
        return Err(HenselLiftError::ZeroTargetLevel);
    }

    let mut modulus = BigUint::from(2u8);
    let mut roots = vec![positive_mod_biguint(value, &modulus)];

    for _level in 2..=e {
        let previous_modulus = modulus.clone();
        modulus <<= 1usize;
        let target = positive_mod_biguint(value, &modulus);
        roots = lift_two_power_roots_one_level(&roots, &previous_modulus, &modulus, &target);

        if roots.is_empty() {
            return Err(HenselLiftError::NoSquareRootModuloPrimePower);
        }
    }

    Ok(roots)
}

fn lift_two_power_roots_one_level(
    roots: &[BigUint],
    previous_modulus: &BigUint,
    modulus: &BigUint,
    target: &BigUint,
) -> Vec<BigUint> {
    let mut lifted = Vec::new();

    for root in roots {
        let candidates = [root.clone(), root + previous_modulus];
        for candidate in candidates {
            if (&candidate * &candidate) % modulus == *target {
                lifted.push(candidate);
            }
        }
    }

    lifted.sort();
    lifted.dedup();
    lifted
}
