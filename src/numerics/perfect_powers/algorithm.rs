use num_bigint::BigUint;

use crate::numerics::perfect_powers::{
    PerfectPowerOutcome, PerfectPowerReport, PerfectPowerSearchConfig, candidate::detect_candidate,
    exponents::prime_exponents_through, validation::validate_search_input,
};

/// Detects whether `N` is a perfect power `a^q` with prime exponent `q`.
///
/// It assumes `N > 1` and `gcd(N, 6) = 1`. For each prime `q ≤ ⌊log₂ N⌋`, it builds
/// `f_q(x) = x^q − N`, chooses a Hensel prime that does not divide `q`, recovers
/// certified integer roots of `f_q`, and verifies `a^q = N` by exact exponentiation.
///
/// The Hensel prime choice is:
///
/// 1. for `q = 2`, use `p = 3`;
/// 2. for odd prime `q`, use `p = 2`.
///
/// Complexity: let `n = ⌈log₂ N⌉` and let `M(n)` be the cost of multiplying
/// `n`-bit integers. The algorithm tests `π(n) = Θ(n/log n)` prime exponents.
/// For one exponent `q`, the sparse polynomial `xᵠ − N` has `O(1)` terms,
/// evaluation costs `O(log q · M(n))`, and the Hensel precision is `Θ(n/q)`.
/// Summing over prime `q ≤ n` gives `O(n log n · M(n))` bit operations, plus
/// lower-order sieve and exact-power-check work. With quasi-linear integer
/// multiplication, this is quasi-quadratic in `n`.
pub(crate) fn detect_perfect_power(
    input: &BigUint,
    config: PerfectPowerSearchConfig,
) -> PerfectPowerReport {
    let max_exponent = match validate_search_input(input) {
        Ok(max_exponent) => max_exponent,
        Err(outcome) => {
            return PerfectPowerReport::new(input.clone(), config, Vec::new(), outcome);
        }
    };

    let mut candidate_reports = Vec::new();
    for exponent in prime_exponents_through(max_exponent) {
        let candidate = match detect_candidate(input, &config, exponent) {
            Ok(candidate) => candidate,
            Err(outcome) => {
                return PerfectPowerReport::new(input.clone(), config, candidate_reports, outcome);
            }
        };

        let (report, outcome) = candidate;
        candidate_reports.push(report);
        if let Some(outcome) = outcome {
            return PerfectPowerReport::new(input.clone(), config, candidate_reports, outcome);
        }
    }

    PerfectPowerReport::new(
        input.clone(),
        config,
        candidate_reports,
        PerfectPowerOutcome::NotPerfectPower,
    )
}
