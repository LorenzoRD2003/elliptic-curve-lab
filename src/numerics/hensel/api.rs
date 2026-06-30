use num_bigint::{BigInt, BigUint};
use num_prime::nt_funcs::is_prime;
use num_traits::Zero;

use super::{
    error::HenselLiftError,
    polynomial::{
        evaluate_derivative, evaluate_polynomial, positive_mod_bigint, positive_mod_biguint,
    },
    step::HenselLiftStep,
    trace::HenselLiftTrace,
};
use crate::numerics::inverse_mod_biguint;

/// Lifts one simple root from modulo `p^k` to modulo `p^(k+1)`.
///
/// The polynomial is represented by coefficients in ascending degree order:
/// `coefficients[i]` is the coefficient of `x^i`.
///
/// If `f(x_k) = p^k a`, the simple-root formula chooses the unique digit
/// `t mod p` satisfying
///
/// `a + t f'(x_k) = 0 mod p`,
///
/// and returns `x_{k+1} = x_k + t p^k`.
///
/// Complexity: `Θ(d)` exact integer multiplications and additions for a
/// degree-`d` polynomial, plus one modular inverse modulo `p`.
pub(crate) fn hensel_lift_simple_root_step(
    coefs: &[BigInt],
    root_mod_p_to_k: &BigInt,
    p: &BigUint,
    level: u32,
) -> Result<HenselLiftStep, HenselLiftError> {
    validate_simple_hensel_input(coefs, p)?;

    if level == 0 {
        return Err(HenselLiftError::ZeroTargetLevel);
    }

    let p_to_k = BigInt::from(p.pow(level));
    let value = evaluate_polynomial(coefs, root_mod_p_to_k);
    if !positive_mod_bigint(&value, &p_to_k).is_zero() {
        return Err(HenselLiftError::RootDoesNotSolveCurrentModulus);
    }

    let derivative = evaluate_derivative(coefs, root_mod_p_to_k);
    let derivative_mod_p = positive_mod_biguint(&derivative, p);
    if derivative_mod_p.is_zero() {
        return Err(HenselLiftError::SingularDerivativeModPrime);
    }

    let inverse = inverse_mod_biguint(&derivative_mod_p, p)
        .ok_or(HenselLiftError::SingularDerivativeModPrime)?;
    let quotient_mod_p = positive_mod_biguint(&(value / &p_to_k), p);
    let lift_digit = if quotient_mod_p.is_zero() {
        BigUint::zero()
    } else {
        (p - ((quotient_mod_p * inverse) % p)) % p
    };
    let root_after = root_mod_p_to_k + BigInt::from(lift_digit.clone()) * p_to_k;

    Ok(HenselLiftStep::new(
        level,
        root_mod_p_to_k.clone(),
        lift_digit,
        root_after,
    ))
}

/// Repeatedly applies simple-root Hensel lifting from modulo `p` to modulo
/// `p^e`.
///
/// The starting root must solve `f(x) = 0 mod p`, and each step records the
/// correction digit in `x_{k+1} = x_k + t p^k`.
///
/// Complexity: `Θ(e * d)` exact polynomial evaluation work for a degree-`d`
/// polynomial, plus one modular inverse modulo `p` per step.
pub(crate) fn hensel_lift_simple_root(
    coefs: &[BigInt],
    root_mod_p: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<HenselLiftTrace, HenselLiftError> {
    validate_simple_hensel_input(coefs, p)?;

    if e == 0 {
        return Err(HenselLiftError::ZeroTargetLevel);
    }

    let initial_root = normalize_root(root_mod_p, p);
    if e == 1 {
        ensure_root_solves_level(coefs, &initial_root, p, 1)?;
        return Ok(HenselLiftTrace::new(
            p.clone(),
            coefs.to_vec(),
            initial_root,
            Vec::new(),
        ));
    }

    let mut root = initial_root.clone();
    let mut steps = Vec::with_capacity(e.saturating_sub(1) as usize);
    for level in 1..e {
        let step = hensel_lift_simple_root_step(coefs, &root, p, level)?;
        root = step.root_after().clone();
        steps.push(step);
    }

    Ok(HenselLiftTrace::new(
        p.clone(),
        coefs.to_vec(),
        initial_root,
        steps,
    ))
}

pub(super) fn validate_simple_hensel_input(
    coefficients: &[BigInt],
    prime: &BigUint,
) -> Result<(), HenselLiftError> {
    if coefficients.is_empty() {
        return Err(HenselLiftError::EmptyPolynomial);
    }
    if coefficients.len() == 1 {
        return Err(HenselLiftError::ConstantPolynomial);
    }
    if prime < &BigUint::from(2u8) || !is_prime(prime, None).probably() {
        return Err(HenselLiftError::NonPrimeModulus);
    }
    Ok(())
}

fn ensure_root_solves_level(
    coefficients: &[BigInt],
    root: &BigInt,
    prime: &BigUint,
    level: u32,
) -> Result<(), HenselLiftError> {
    let modulus = BigInt::from(prime.pow(level));
    let value = evaluate_polynomial(coefficients, root);
    if positive_mod_bigint(&value, &modulus).is_zero() {
        Ok(())
    } else {
        Err(HenselLiftError::RootDoesNotSolveCurrentModulus)
    }
}

pub(super) fn normalize_root(root: &BigInt, prime: &BigUint) -> BigInt {
    positive_mod_bigint(root, &BigInt::from(prime.clone()))
}
