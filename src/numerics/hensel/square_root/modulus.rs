use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

use crate::numerics::{
    chinese_remainder::{Congruence, solve_coprime_congruences},
    hensel::HenselLiftError,
    positive_mod_biguint,
    prime_powers::NormalizedPrimePowerFactorization,
};

use super::{
    sqrt_mod_odd_prime_power, sqrt_mod_odd_prime_power_divisible_radicand, sqrt_mod_two_power,
};

/// Returns all square roots of `value` modulo an arbitrary integer `m ≥ 2`.
///
/// The route factors `m = Πᵢ pᵢ^eᵢ`, solves the square-root problem modulo each
/// prime power `pᵢ^eᵢ`, and combines every compatible tuple of local roots with
/// CRT. The returned roots are sorted canonical representatives modulo `m`.
///
/// Complexity: the cost of factoring `m`, plus the sum of local prime-power
/// square-root costs, plus `Θ(R · r²)` gcd checks and `Θ(R · r)` CRT
/// combinations, where `r` is the number of distinct prime factors of `m` and
/// `R = Πᵢ |Roots(a mod pᵢ^eᵢ)|` is the number of returned roots.
pub(crate) fn sqrt_mod_m(value: &BigInt, m: &BigUint) -> Result<Vec<BigUint>, HenselLiftError> {
    if m < &BigUint::from(2u8) {
        return Err(HenselLiftError::TrivialModulus);
    }

    let factors = NormalizedPrimePowerFactorization::factor(m)
        .map_err(|_| HenselLiftError::TrivialModulus)?
        .into_factors();
    let local_roots = factors
        .iter()
        .map(|(p, e)| square_roots_mod_prime_power(value, p, *e))
        .collect::<Result<Vec<_>, _>>()?;

    combine_local_square_roots(&local_roots)
}

fn square_roots_mod_prime_power(
    value: &BigInt,
    p: &BigUint,
    e: u32,
) -> Result<(BigUint, Vec<BigUint>), HenselLiftError> {
    let prime_power = p.pow(e);
    let roots = if p == &BigUint::from(2u8) {
        sqrt_mod_two_power(value, e)?
    } else if positive_mod_biguint(value, p).is_zero() {
        sqrt_mod_odd_prime_power_divisible_radicand(value, p, e)?
    } else {
        let (left, right) = sqrt_mod_odd_prime_power(value, p, e).map_err(|error| match error {
            HenselLiftError::QuadraticNonResidueModPrime => {
                HenselLiftError::NoSquareRootModuloPrimePower
            }
            other => other,
        })?;
        vec![left, right]
    };

    Ok((prime_power, roots))
}

fn combine_local_square_roots(
    local_roots: &[(BigUint, Vec<BigUint>)],
) -> Result<Vec<BigUint>, HenselLiftError> {
    let mut roots = Vec::new();
    let mut congruences = Vec::with_capacity(local_roots.len());
    collect_crt_combinations(local_roots, 0, &mut congruences, &mut roots)?;
    roots.sort();
    roots.dedup();
    Ok(roots)
}

fn collect_crt_combinations(
    local_roots: &[(BigUint, Vec<BigUint>)],
    index: usize,
    congruences: &mut Vec<Congruence>,
    roots: &mut Vec<BigUint>,
) -> Result<(), HenselLiftError> {
    if index == local_roots.len() {
        let solution = solve_coprime_congruences(congruences)
            .expect("distinct prime powers should be pairwise coprime");
        roots.push(solution.residue().clone());
        return Ok(());
    }

    let (modulus, residues) = &local_roots[index];
    for residue in residues {
        congruences.push(
            Congruence::new(residue.clone(), modulus.clone())
                .expect("prime powers should be valid CRT moduli"),
        );
        collect_crt_combinations(local_roots, index + 1, congruences, roots)?;
        congruences.pop();
    }

    Ok(())
}
