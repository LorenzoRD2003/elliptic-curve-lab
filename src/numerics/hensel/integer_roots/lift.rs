use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

use crate::{
    numerics::hensel::{
        HenselIntegerRootTrace, HenselLiftError, api::hensel_lift_simple_root_for_polynomial,
        polynomial::positive_mod_bigint,
    },
    polynomials::IntegerPolynomial,
};

/// Lifts a simple root modulo `p` and certifies the corresponding integer root.
///
/// The polynomial is represented by an [`IntegerPolynomial`].
/// The supplied `root_bound` is a promised absolute bound `B₀` for the target
/// integer root. The helper chooses the least exponent `e` with `pᵉ > 2B₀`,
/// lifts the seed to modulo `pᵉ`, centers the residue in the symmetric
/// interval around zero, and finally verifies `f(r) = 0` over `ℤ`.
///
/// This first integer-root surface intentionally handles only the simple-root
/// Hensel case where `f′(x₀)` is a unit modulo `p`.
///
/// Correctness intuition: a class modulo `pᵉ` can contain at most one integer in
/// `[-B₀, B₀]` once `pᵉ > 2B₀`. The Hensel lift supplies the class; the final
/// exact evaluation supplies the certificate that the centered representative
/// is a genuine root in `ℤ`.
///
/// Complexity: let `s` be the number of non-zero terms, `d = deg f`, and
/// `e = ⌈logₚ(2B₀ + 1)⌉`. This route performs `e - 1` simple Hensel steps plus
/// one final exact evaluation. With the sparse evaluator, that is
/// `Θ(e·s·log d)` exact integer-multiplication work on integers whose size grows
/// up to `O(log B₀)` bits, plus one modular inverse modulo `p` per Hensel step.
pub(crate) fn hensel_lift_integer_root(
    polynomial: &IntegerPolynomial,
    root_mod_p: &BigInt,
    prime: &BigUint,
    root_bound: &BigUint,
) -> Result<HenselIntegerRootTrace, HenselLiftError> {
    let target_level = target_level_for_root_bound(prime, root_bound)?;
    let lift_trace =
        hensel_lift_simple_root_for_polynomial(polynomial, root_mod_p, prime, target_level)?;
    let modulus = prime.pow(target_level);
    let candidate_root = centered_representative(lift_trace.final_root(), &modulus);

    if candidate_root.magnitude() > root_bound || !polynomial.evaluate(&candidate_root).is_zero() {
        return Err(HenselLiftError::IntegerRootNotCertifiedInBound);
    }

    Ok(HenselIntegerRootTrace::new(
        lift_trace,
        root_bound.clone(),
        modulus,
        candidate_root,
    ))
}

fn target_level_for_root_bound(
    prime: &BigUint,
    root_bound: &BigUint,
) -> Result<u32, HenselLiftError> {
    let target = root_bound << 1usize;
    let mut level = 1u32;
    let mut modulus = prime.clone();

    while modulus <= target {
        level = level
            .checked_add(1)
            .ok_or(HenselLiftError::TargetLevelOverflow)?;
        modulus *= prime;
    }
    Ok(level)
}

fn centered_representative(root: &BigInt, modulus: &BigUint) -> BigInt {
    let modulus_bigint = BigInt::from(modulus.clone());
    let residue = positive_mod_bigint(root, &modulus_bigint);
    let doubled = &residue << 1usize;

    if doubled > modulus_bigint {
        residue - modulus_bigint
    } else {
        residue
    }
}
