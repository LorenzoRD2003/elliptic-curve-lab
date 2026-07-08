use num_bigint::{BigInt, BigUint};
use num_prime::nt_funcs::is_prime;
use num_traits::{ToPrimitive, Zero};

use crate::{
    numerics::{
        hensel::{
            HenselIntegerRootSearchConfig, HenselIntegerRootSearchReport, HenselLiftError,
            hensel_lift_integer_root,
        },
        positive_mod_biguint,
    },
    polynomials::IntegerPolynomial,
};

/// Finds all bounded integer roots certified from simple roots modulo one prime.
///
/// The algorithm scans every residue `x₀ mod p`. If `f(x₀) ≠ 0 mod p`, the
/// residue cannot be a seed and is ignored. If `f(x₀) = 0 mod p` but
/// `f′(x₀) = 0 mod p`, the residue is counted as singular and skipped by this
/// first simple-root route. Every remaining seed is lifted with
/// [`hensel_lift_integer_root`].
///
/// Why this can certify roots:
///
/// 1. Hensel lifting gives a unique lift modulo `pᵉ` for each simple seed.
/// 2. The chosen precision satisfies `pᵉ > 2B₀`, so a residue class modulo
///    `pᵉ` contains at most one integer in `[-B₀, B₀]`.
/// 3. The centered candidate is checked by exact evaluation `f(r) = 0` in
///    `ℤ`, so a p-adic lift is never mistaken for an integer root.
///
/// Therefore an empty result is a complete negative result for simple seeds
/// modulo this prime. It does not yet rule out roots that are singular modulo
/// this prime; callers can choose another prime or wait for a later singular
/// lifting route.
///
/// Complexity: the seed scan is linear in the prime `p`: it evaluates `f` and
/// `f′` at each residue modulo `p`, so it costs `Θ(p·d)` field-sized work for
/// `d = deg f`. Each simple seed that survives the scan is then lifted and
/// certified by [`hensel_lift_integer_root`].
pub(crate) fn find_integer_roots_by_hensel(
    polynomial: &IntegerPolynomial,
    config: HenselIntegerRootSearchConfig,
) -> Result<HenselIntegerRootSearchReport, HenselLiftError> {
    validate_integer_root_search_input(polynomial, &config)?;

    let prime_u64 = config
        .prime()
        .to_u64()
        .ok_or(HenselLiftError::SeedScanLimitExceeded)?;

    if prime_u64 > config.max_seed_scan() {
        return Err(HenselLiftError::SeedScanLimitExceeded);
    }

    let mut simple_seed_count = 0usize;
    let mut singular_seed_count = 0usize;
    let mut uncertified_seed_count = 0usize;
    let mut traces = Vec::new();
    let mut certified_roots = Vec::new();

    for seed in 0..prime_u64 {
        let seed = BigInt::from(seed);
        if !positive_mod_biguint(&polynomial.evaluate(&seed), config.prime()).is_zero() {
            continue;
        }

        if positive_mod_biguint(&polynomial.evaluate_derivative(&seed), config.prime()).is_zero() {
            singular_seed_count += 1;
            continue;
        }

        simple_seed_count += 1;
        match hensel_lift_integer_root(polynomial, &seed, config.prime(), config.root_bound()) {
            Ok(trace) => {
                certified_roots.push(trace.candidate_root().clone());
                traces.push(trace);
            }
            Err(HenselLiftError::IntegerRootNotCertifiedInBound) => {
                uncertified_seed_count += 1;
            }
            Err(error) => return Err(error),
        }
    }

    certified_roots.sort();
    certified_roots.dedup();

    Ok(HenselIntegerRootSearchReport::new(
        config,
        simple_seed_count,
        singular_seed_count,
        uncertified_seed_count,
        certified_roots,
        traces,
    ))
}

fn validate_integer_root_search_input(
    polynomial: &IntegerPolynomial,
    config: &HenselIntegerRootSearchConfig,
) -> Result<(), HenselLiftError> {
    if polynomial.is_zero() {
        return Err(HenselLiftError::EmptyPolynomial);
    }
    if polynomial.is_constant() {
        return Err(HenselLiftError::ConstantPolynomial);
    }
    if config.prime() < &BigUint::from(2u8) || !is_prime(config.prime(), None).probably() {
        return Err(HenselLiftError::NonPrimeModulus);
    }
    Ok(())
}
