use num_bigint::{BigInt, BigUint};

use crate::numerics::{
    hensel::{HenselLiftError, hensel_lift_square_root_fast},
    positive_mod_biguint,
};

use super::{tonelli::tonelli_shanks_mod_odd_prime, validation::validate_unit_input};

/// Returns the two square roots of `value` modulo `p^e` for odd prime `p`.
///
/// - `p` must be an odd prime
/// - `e ≥ 1`
/// - `value` must be non-zero modulo `p`
/// - `value` must be a quadratic residue modulo `p`
///
/// The algorithm reduces `value` modulo `p`, finds a root modulo `p` with
/// Tonelli-Shanks, lifts it to `p^e` with the fast square-root Hensel route,
/// and returns `(r, -r mod p^e)` as canonical representatives.
///
/// Complexity: dominated by Tonelli-Shanks modulo `p` plus fast Hensel lifting
/// to `p^e`.
pub(crate) fn sqrt_mod_odd_prime_power(
    value: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<(BigUint, BigUint), HenselLiftError> {
    validate_unit_input(value, p, e)?;

    let root_mod_p = tonelli_shanks_mod_odd_prime(&positive_mod_biguint(value, p), p)
        .ok_or(HenselLiftError::QuadraticNonResidueModPrime)?;
    let trace = hensel_lift_square_root_fast(value, &BigInt::from(root_mod_p), p, e)?;
    let modulus = p.pow(e);
    let root = positive_mod_biguint(trace.final_root(), &modulus);
    let neg_root = (&modulus - &root) % &modulus;

    if root <= neg_root {
        Ok((root, neg_root))
    } else {
        Ok((neg_root, root))
    }
}
