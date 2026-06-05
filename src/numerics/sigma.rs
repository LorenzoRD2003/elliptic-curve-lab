use std::num::NonZeroUsize;

use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Returns the classical divisor-power sum `σ_power(n) = Σ_{d | n} d^power`
/// by scanning trial divisors up to `⌊√n⌋`.
///
/// This is the most literal educational implementation in the crate:
/// whenever `d | n`, it adds both `d^power` and `(n / d)^power`, taking care
/// not to double-count the square-root divisor when `n` is a square.
///
/// Complexity:
/// - `Θ(√n)` divisibility checks
/// - up to two `BigInt` exponentiations for each divisor pair that is found
/// - `Θ(1)` auxiliary memory beyond the returned integer
pub fn sigma_power_sum_naive(n: NonZeroUsize, power: u32) -> BigInt {
    let mut sum = BigInt::zero();
    let mut divisor = 1usize;
    let n = n.get();
    while divisor * divisor <= n {
        if n.is_multiple_of(divisor) {
            let paired_divisor = n / divisor;
            sum += BigInt::from(divisor).pow(power);

            if paired_divisor != divisor {
                sum += BigInt::from(paired_divisor).pow(power);
            }
        }
        divisor += 1;
    }
    sum
}

/// Returns the classical divisor-power sum
/// `σ_power(n) = Σ_{d | n} d^power`
/// by first factoring `n` and then using multiplicativity.
///
/// If `n = ∏ p_i^{e_i}`, then
/// `σ_power(n) = ∏ (1 + p_i^power + p_i^{2 power} + ... + p_i^{e_i power})`.
///
/// This avoids summing over all divisors explicitly. Once the prime-power
/// factorization is known, the remaining work is only the product of short
/// geometric sums attached to the prime powers.
///
/// In this first educational implementation the factorization itself is still
/// obtained by trial division, so the dominant cost remains:
///
/// - `Θ(√n)` trial divisions in the worst case
/// - at most `O(log n)` geometric-series steps after factoring, since the
///   total prime-factor multiplicity of `n` is bounded above by `⌊log₂ n⌋`.
pub fn sigma_power_sum_factorized(n: NonZeroUsize, power: u32) -> BigInt {
    trial_division_prime_power_factorization(n)
        .into_iter()
        .map(|factor| sigma_prime_power_sum(factor.prime, factor.exponent, power))
        .product()
}

/// Returns the whole table `σ_power(0), σ_power(1), ..., σ_power(bound)`
/// using a divisor-sieve traversal.
///
/// The classical divisor sum is only defined for positive integers, so the
/// first entry is stored as the explicit sentinel value `0`. For every
/// positive divisor `d ≤ bound`, the algorithm computes `d^power` once and
/// adds it to all multiples of `d`.
///
/// Mathematically this is the identity `σ_power(m) = Σ_{d | m} d^power`
/// turned inside out: instead of fixing `m` and enumerating its divisors, we
/// fix a divisor `d` and contribute it to every multiple that contains it.
///
/// - `Θ(bound log bound)` `BigInt` additions, since
///   `Σ_{d=1}^{bound} ⌊bound / d⌋ = Θ(bound ln bound)`
/// - one `BigInt` exponentiation per divisor `d`
/// - `Θ(bound)` auxiliary memory for the output table
///
/// This is the preferred implementation when a caller needs many consecutive
/// sigma values, for example when assembling truncated Eisenstein-series
/// coefficients up to some `q^N`.
pub fn sigma_power_sums_up_to(bound: usize, power: u32) -> Vec<BigInt> {
    let mut sums = vec![BigInt::zero(); bound + 1];

    for divisor in 1..=bound {
        let divisor_power = BigInt::from(divisor).pow(power);
        let mut multiple = divisor;

        while multiple <= bound {
            sums[multiple] += &divisor_power;
            multiple += divisor;
        }
    }

    sums
}

/// A single prime-power term `p^e` appearing in the factorization of `n`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PrimePowerFactor {
    prime: usize,
    exponent: u32,
}

/// Factors a positive integer into prime powers by trial division.
///
/// If `n = ∏ p_i^{e_i}`, this returns the list of pairs `(p_i, e_i)`
/// in increasing prime order.
///
/// - `Θ(√n)` trial divisions in the worst case
/// - at most `O(log n)` stored factors, since the number of distinct prime
///   divisors is bounded above by `⌊log₂ n⌋.
///
/// The implementation checks `2` first and then only odd candidates, which
/// preserves the same asymptotic bound while avoiding obviously redundant
/// even-divisor tests.
fn trial_division_prime_power_factorization(n: NonZeroUsize) -> Vec<PrimePowerFactor> {
    let mut remaining = n.get();
    let mut factors = Vec::new();
    let mut divisor = 2usize;

    while divisor * divisor <= remaining {
        let mut exponent = 0u32;
        while remaining.is_multiple_of(divisor) {
            remaining /= divisor;
            exponent += 1;
        }

        if exponent > 0 {
            factors.push(PrimePowerFactor {
                prime: divisor,
                exponent,
            });
        }

        divisor = next_trial_divisor(divisor);
    }

    if remaining > 1 {
        factors.push(PrimePowerFactor {
            prime: remaining,
            exponent: 1,
        });
    }

    factors
}

/// Returns the local geometric factor
/// `1 + p^power + p^{2 power} + ... + p^{exponent * power}`.
///
/// This is the prime-power contribution
/// `σ_power(p^exponent)` used in the multiplicative formula for `σ_power(n)`.
///
/// - `Θ(exponent)` `BigInt` multiplications and additions
/// - `Θ(1)` auxiliary memory beyond the returned integer
fn sigma_prime_power_sum(prime: usize, exponent: u32, power: u32) -> BigInt {
    let base = BigInt::from(prime).pow(power);
    let mut sum = BigInt::one();
    let mut current_term = BigInt::one();
    for _ in 0..exponent {
        current_term *= &base;
        sum += &current_term;
    }
    sum
}

/// Returns the next trial divisor after the current one.
///
/// The sequence is `2, 3, 5, 7, 9, ...`: after handling `2` separately, the
/// factorization routine only visits odd candidates.
fn next_trial_divisor(divisor: usize) -> usize {
    if divisor == 2 { 3 } else { divisor + 2 }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use num_bigint::BigInt;

    use super::{
        sigma_power_sum_factorized, sigma_power_sum_naive, sigma_power_sums_up_to,
        sigma_prime_power_sum, trial_division_prime_power_factorization,
    };

    fn positive(n: usize) -> NonZeroUsize {
        NonZeroUsize::new(n).unwrap()
    }

    #[test]
    fn naive_sigma_matches_small_known_values() {
        assert_eq!(sigma_power_sum_naive(positive(1), 0), BigInt::from(1usize));
        assert_eq!(sigma_power_sum_naive(positive(6), 1), BigInt::from(12usize));
        assert_eq!(
            sigma_power_sum_naive(positive(12), 1),
            BigInt::from(28usize)
        );
        assert_eq!(
            sigma_power_sum_naive(positive(12), 2),
            BigInt::from(210usize)
        );
    }

    #[test]
    fn factorization_matches_small_known_values() {
        assert_eq!(
            sigma_power_sum_factorized(positive(6), 1),
            BigInt::from(12usize)
        );
        assert_eq!(
            sigma_power_sum_factorized(positive(12), 1),
            BigInt::from(28usize)
        );
        assert_eq!(
            sigma_power_sum_factorized(positive(12), 2),
            BigInt::from(210usize)
        );
    }

    #[test]
    fn factorization_and_naive_paths_agree_on_small_grid() {
        for n in 1..=32 {
            for power in 0..=6 {
                assert_eq!(
                    sigma_power_sum_factorized(positive(n), power),
                    sigma_power_sum_naive(positive(n), power)
                );
            }
        }
    }

    #[test]
    fn divisor_sieve_matches_naive_values() {
        let table = sigma_power_sums_up_to(20, 3);

        assert_eq!(table[0], BigInt::from(0usize));
        for (n, sigma_value) in table.iter().enumerate().take(21).skip(1) {
            assert_eq!(*sigma_value, sigma_power_sum_naive(positive(n), 3));
        }
    }

    #[test]
    fn trial_division_factorization_keeps_prime_powers() {
        assert_eq!(
            trial_division_prime_power_factorization(positive(360)),
            vec![
                super::PrimePowerFactor {
                    prime: 2,
                    exponent: 3,
                },
                super::PrimePowerFactor {
                    prime: 3,
                    exponent: 2,
                },
                super::PrimePowerFactor {
                    prime: 5,
                    exponent: 1,
                },
            ]
        );
    }

    #[test]
    fn prime_power_sigma_component_matches_geometric_sum() {
        assert_eq!(sigma_prime_power_sum(2, 3, 1), BigInt::from(15usize));
        assert_eq!(sigma_prime_power_sum(3, 2, 2), BigInt::from(91usize));
    }

    #[test]
    fn divisor_sieve_handles_zero_bound_with_a_zero_sentinel() {
        assert_eq!(sigma_power_sums_up_to(0, 5), vec![BigInt::from(0usize)]);
    }
}
