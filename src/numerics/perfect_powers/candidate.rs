use num_bigint::{BigInt, BigUint, Sign};
use num_traits::One;

use crate::numerics::{
    hensel::{HenselIntegerRootSearchConfig, find_integer_roots_by_hensel},
    perfect_powers::{
        PerfectPowerCandidateReport, PerfectPowerOutcome, PerfectPowerSearchConfig,
        exponents::hensel_prime_for_exponent,
    },
};
use crate::polynomials::IntegerPolynomial;

/// Tests one prime exponent `q` in the perfect-power search.
///
/// For fixed `N` and `q`, this helper turns `N = a^q` into the integer-root
/// problem `f_q(x) = x^q − N`. It then:
///
/// 1. builds the integer polynomial `f_q`;
/// 2. computes a safe root bound `B = 2^⌈bits(N)/q⌉`;
/// 3. chooses the staged Hensel prime for `q`;
/// 4. asks the bounded integer-root search to certify roots of `f_q`;
/// 5. accepts a positive certified root only after the exact check `aᵠ = N`.
///
/// The returned candidate report records the Hensel prime, the root bound, and
/// every certified integer root. The optional outcome is
/// `Some(PerfectPower { base, exponent: q })` exactly when this exponent certifies
/// a decomposition of `N`; otherwise it is `None`.
///
/// Complexity: building sparse `xᵠ − N` costs `Θ(1)`. The dominant work is
/// delegated to [`find_integer_roots_by_hensel`]. After that, this helper
/// performs one exact exponentiation check for each certified positive root.
pub(super) fn detect_candidate(
    input: &BigUint,
    config: &PerfectPowerSearchConfig,
    exponent: u32,
) -> Result<(PerfectPowerCandidateReport, Option<PerfectPowerOutcome>), PerfectPowerOutcome> {
    let polynomial = IntegerPolynomial::x_power_minus_constant(input, exponent)
        .ok_or(PerfectPowerOutcome::ExponentTooLarge)?;
    let root_bound = root_bound_for_exponent(input, u64::from(exponent))
        .ok_or(PerfectPowerOutcome::ExponentTooLarge)?;
    let hensel_prime = hensel_prime_for_exponent(exponent);
    let hensel_config =
        HenselIntegerRootSearchConfig::new(hensel_prime.clone(), root_bound.clone())
            .with_max_seed_scan(config.max_seed_scan());
    let root_report = find_integer_roots_by_hensel(&polynomial, hensel_config)
        .map_err(PerfectPowerOutcome::HenselFailure)?;

    let certified_roots = root_report.certified_roots().to_vec();
    let outcome = certified_roots
        .iter()
        .filter_map(positive_base_from_root)
        .find(|base| base.pow(exponent) == *input)
        .map(|base| PerfectPowerOutcome::PerfectPower { base, exponent });

    Ok((
        PerfectPowerCandidateReport::new(exponent, hensel_prime, root_bound, certified_roots),
        outcome,
    ))
}

fn positive_base_from_root(root: &BigInt) -> Option<BigUint> {
    match root.sign() {
        Sign::Minus | Sign::NoSign => None,
        Sign::Plus => root.to_biguint(),
    }
}

/// Returns a safe root bound for `x^q − N`.
///
/// If `N` has `b` bits, then `N < 2^b`. Any positive integer solution `a^q = N`
/// therefore satisfies `a < 2^(b/q)`. This helper returns the slightly
/// generous bound `B = 2^⌈b/q⌉`, which is enough for the downstream integer-root
/// search to certify or reject roots inside `|x| ≤ B`.
///
/// Complexity: `Θ(1)` arithmetic, plus construction of a `BigUint` with `Θ(⌈b/q⌉)` bits.
fn root_bound_for_exponent(input: &BigUint, exponent: u64) -> Option<BigUint> {
    let bits = input.bits();
    let bound_bits = bits.div_ceil(exponent);
    let shift = usize::try_from(bound_bits).ok()?;
    Some(BigUint::one() << shift)
}
