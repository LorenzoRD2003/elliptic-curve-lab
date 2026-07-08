use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

use crate::numerics::hensel::{
    HenselLiftError, HenselLiftTrace, HenselSquareRootFastStep, HenselSquareRootFastTrace,
    api::{normalize_root, validate_simple_hensel_input},
    hensel_lift_simple_root,
};
use crate::numerics::{inverse_mod_biguint, positive_mod_bigint, positive_mod_biguint};

/// Lifts one simple square root from modulo `p` to modulo `p^e`.
///
/// This is the specialized route for the equation `x^2 = value`. It delegates
/// to the general simple-root Hensel engine with the polynomial
/// `f(x) = x^2 - value`, so it still requires the simple-root condition
/// `2x != 0 mod p`.
///
/// Complexity: `Θ(e)` exact integer work, plus one modular inverse modulo
/// `p` per lifting step.
pub(crate) fn hensel_lift_square_root(
    value: &BigInt,
    root_mod_p: &BigInt,
    prime: &BigUint,
    e: u32,
) -> Result<HenselLiftTrace, HenselLiftError> {
    hensel_lift_simple_root(&square_root_polynomial(value), root_mod_p, prime, e)
}

/// Lifts one simple square root by Newton-Hensel precision doubling:
///
/// `p -> p^2 -> p^4 -> p^8 -> ... -> p^e`.
///
/// For a current solution modulo `q = p^k`, the next representative is
///
/// `x_new = x_old - (x_old^2 - value) * (2x_old)^(-1) mod q^2`.
///
/// The final step uses `p^e` directly when the requested level lies between
/// two powers reached by doubling. This route is intentionally specialized to
/// square roots and requires `2x_old` to be a unit modulo `p`.
///
/// Complexity: `Θ(log e)` Newton-Hensel steps. Step `i` performs Θ(1)
/// operations modulo roughly `p^(2^i)`, including one modular inverse. The total
/// bit-complexity is dominated by the final inverse modulo `p^e`.
pub(crate) fn hensel_lift_square_root_fast(
    value: &BigInt,
    root_mod_p: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<HenselSquareRootFastTrace, HenselLiftError> {
    validate_simple_hensel_input(&square_root_polynomial(value), p)?;

    if e == 0 {
        return Err(HenselLiftError::ZeroTargetLevel);
    }

    let initial_root = normalize_root(root_mod_p, p);
    let prime_modulus = BigInt::from(p.clone());
    if !positive_mod_bigint(&(&initial_root * &initial_root - value), &prime_modulus).is_zero() {
        return Err(HenselLiftError::RootDoesNotSolveCurrentModulus);
    }

    if e == 1 {
        return Ok(HenselSquareRootFastTrace::new(
            p.clone(),
            value.clone(),
            initial_root,
            e,
            Vec::new(),
        ));
    }

    let mut current_level = 1;
    let mut root = initial_root.clone();
    let mut steps = Vec::new();

    while current_level < e {
        let next_level = current_level.saturating_mul(2).min(e);
        let root_before = root.clone();
        root = fast_square_root_step(value, &root, p, next_level)?;
        steps.push(HenselSquareRootFastStep::new(
            current_level,
            next_level,
            root_before,
            root.clone(),
        ));
        current_level = next_level;
    }

    Ok(HenselSquareRootFastTrace::new(
        p.clone(),
        value.clone(),
        initial_root,
        e,
        steps,
    ))
}

pub(super) fn square_root_polynomial(value: &BigInt) -> [BigInt; 3] {
    [-value, BigInt::zero(), BigInt::from(1u8)]
}

fn fast_square_root_step(
    value: &BigInt,
    root: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<BigInt, HenselLiftError> {
    let modulus = p.pow(e);
    let modulus_bigint = BigInt::from(modulus.clone());
    let residual = root * root - value;
    let derivative = BigInt::from(2u8) * root;
    let derivative_modulus = positive_mod_biguint(&derivative, &modulus);
    let inverse = inverse_mod_biguint(&derivative_modulus, &modulus)
        .ok_or(HenselLiftError::SingularDerivativeModPrime)?;
    let correction = positive_mod_bigint(&(residual * BigInt::from(inverse)), &modulus_bigint);

    Ok(positive_mod_bigint(&(root - correction), &modulus_bigint))
}
